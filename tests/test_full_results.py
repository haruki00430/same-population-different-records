"""Run after simulations; skipped only when saved full outputs are absent."""
import json
import numpy as np
import pandas as pd
import pytest
from src.common import ROOT,config
from src.risk_adjustment import recorded_expected

@pytest.fixture
def results():
    if not (ROOT/'results/main/risk_paired.csv').exists():pytest.skip('Full run not executed yet')
    return pd.read_csv(ROOT/'results/main/ac.csv'),pd.read_csv(ROOT/'results/main/risk_paired.csv')

def test_full_analytic_agreement(results):
    ac,_=results
    assert ac.analytic_z_delta.abs().max()<6
    assert (ac.R==10000).all()

def test_risk_expected_relative_response(results):
    _,r=results
    # Population ratio differs slightly from finite-sample E(O/E_B/O/E_A).
    # 0.3 pp is a prespecified numerical validation tolerance, not empirical calibration.
    assert np.max(np.abs(r.mean_RC-r.expected_RC))<.003
    equal=r[np.isclose(r.c_B,.7)&np.isclose(r.g,0)]
    assert np.all(equal.delta==0)
    assert (r.R==10000).all()

def test_risk_development_model_frozen():
    path=ROOT/'results/model_parameters.json'
    if not path.exists():pytest.skip('Model not fit yet')
    model=json.loads(path.read_text());assert model['development_N']==250000
    assert model['gradient_inf']<1e-7
    assert model['B1']!=model['B2']

def test_expected_B1_MC_independent_integration():
    from scipy.special import expit
    cfg=config();path=ROOT/'results/model_parameters.json'
    if not path.exists():pytest.skip('Model not fit yet')
    model=json.loads(path.read_text());rg=np.random.default_rng(920031)
    n=400000;x=rg.normal(size=(n,2));c=rg.random((n,3))<cfg['risk']['comorbidity_prevalence']
    for name in ['B1','B2']:
        co=np.array(model[name]);coded=rg.random(c.shape)<np.where(c,.7,.005)
        p=expit(co[0]+x@co[1:3]+coded@co[3:]);expected=recorded_expected(co,.7,.995,cfg)
        assert abs(p.mean()-expected)<6*p.std(ddof=1)/np.sqrt(n)

def test_spc_output_and_roots():
    if not (ROOT/'tables/table2_tipping_points.csv').exists():pytest.skip('Outputs not built yet')
    t=pd.read_csv(ROOT/'tables/table2_tipping_points.csv');assert t.root_residual.abs().max()<1e-9
    s=pd.read_csv(ROOT/'results/main/spc.csv');assert (s.R==10000).all()
    assert (s.detected_R+s.nondetected_R==s.R).all()
    assert s.signal_probability.between(0,1).all()
    assert (s.wilson_lo<=s.signal_probability).all() and (s.wilson_hi>=s.signal_probability).all()

def test_larger_N_precision_and_persistent_bias():
    path=ROOT/'tables/S8_precision.csv'
    if not path.exists():pytest.skip('Outputs not built yet')
    d=pd.read_csv(path).sort_values('N')
    assert d.delta_sd.iloc[-1]<d.delta_sd.iloc[0]/3
    assert np.all(np.abs(d.delta-.005)<6*d.delta_mcse)
    rel=pd.read_csv(ROOT/'results/sensitivity/reliability.csv')
    for _,group in rel.groupby('record'):
        assert np.max(np.abs(group.estimated_reliability-group.analytic_reliability))<.03
