import numpy as np
from src.common import config,rng_for
from src.dgp import patients
from src.recording import record_patient_view,process_joint
from src.risk_adjustment import calibrate,expected_risk
from src.simulations import safety_cell,process_cell
from src.spc import evaluate_chart

def test_latent_identity_and_immutability():
    cfg=config();rg=rng_for(1,'identity');p=patients(rg,(2000,),calibrate(cfg),cfg)
    u={k:rg.random(shape) for k,shape in [('y',(2000,)),('c',(2000,3)),('e',(2000,)),('a',(2000,))]}
    before=[x.copy() for x in [p.c,p.outcome,p.eligibility,p.action,p.risk]]
    a=record_patient_view(p,.8,.7,.995,.95,.005,.8,u)
    b=record_patient_view(p,1.,1.,.995,.8,.03,1.,u)
    assert a['latent'] is b['latent']
    for x,y in zip(before,[p.c,p.outcome,p.eligibility,p.action,p.risk]):
        np.testing.assert_array_equal(x,y);assert not y.flags.writeable
    assert np.any(a['outcome']!=b['outcome'])
    assert np.all(p.action<=p.eligibility)

def test_analytical_A_and_specificity():
    cfg=config()
    for sp in [1,.995,.99]:
        r=safety_cell(cfg,2000,.05,.8,.9,sp_b=sp)
        assert abs(r['analytic_z_delta'])<6

def test_analytical_C():
    cfg=config()
    for d,f,n in [(.95,.005,.9),(.7,.03,.6),(1.,0,1.)]:
        r=process_cell(cfg,2000,n_b=n,d_b=d,f_b=f)
        assert abs(r['analytic_z_delta'])<6

def test_identity_null_and_seed():
    cfg=config()
    a=safety_cell(cfg,2000,.05,.8,.8)
    assert a['delta']==0 and a['delta_sd']==0
    c=process_cell(cfg,2000)
    assert c['delta']==0 and c['delta_sd']==0
    assert safety_cell(cfg,2000,.05,.8,.9)==safety_cell(cfg,2000,.05,.8,.9)

def test_process_negative_control():
    from src.measures import process_expectation
    assert np.isclose(process_expectation(.3,.75,.7,0,.8),process_expectation(.3,.75,1,0,.8))
    for g in [0,.1,.2]:
        p=process_joint(.3,.75,.75*(1+g),.95,.9,.005,.01,.8,.9)
        assert np.all(p>=-1e-14) and np.isclose(p.sum(),1)

def test_risk_calibration_and_quadrature():
    cfg=config();r=cfg['risk'];a=calibrate(cfg)
    assert abs(expected_risk(a,r['beta'],r['comorbidity_prevalence'])-.08)<1e-10
    assert abs(expected_risk(a,r['beta'],r['comorbidity_prevalence'],100)-.08)<1e-10

def test_spc_rules_and_reset():
    base=np.repeat(.04,24)
    rates=np.array([np.r_[base,np.repeat(.041,24)],np.r_[base,np.repeat(.04,24)],np.r_[base,.07,np.repeat(.04,23)]])
    _,_,_,fr,fs=evaluate_chart(rates,24)
    assert fr[0]==32 and np.isnan(fr[1]) and np.isnan(fr[2])
    assert np.isnan(fs[0]) and fs[2]==25

def test_safety_tipping_identity():
    from src.measures import safety_expectation
    for g in [.1,.2]:
        tipping=.8/(1-g)
        assert np.isclose(safety_expectation(.05,.8),safety_expectation(.05*(1-g),tipping))

def test_conditional_process_expectation_patient_simulation():
    # Independent implementation at patient level checks collapsed category algorithm.
    cfg=config();rg=rng_for(829,'individual');n=500000
    e=rg.random(n)<.3;a=e&(rg.random(n)<.75)
    d=rg.random(n)<np.where(e,.9,.02);num=a&(rg.random(n)<.8)
    value=np.sum(d&num)/np.sum(d)
    from src.measures import process_expectation
    assert abs(value-process_expectation(.3,.75,.9,.02,.8))<.006
