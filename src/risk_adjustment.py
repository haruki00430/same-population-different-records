import itertools
import numpy as np
from scipy.special import expit
from scipy.optimize import brentq, minimize
from numpy.polynomial.hermite import hermgauss
from .common import rng_for, write_json
from .dgp import patients
from .recording import record_binary

def expected_risk(alpha,beta,prevalences,order=60):
    nodes,weights=hermgauss(order)
    normal=nodes*np.sqrt(2*np.sum(np.asarray(beta[:2])**2))
    result=0.
    for bits in itertools.product([0,1],repeat=3):
        prob=np.prod([p if b else 1-p for b,p in zip(bits,prevalences)])
        result+=prob*np.dot(weights,expit(alpha+normal+np.dot(bits,beta[2:])))/np.sqrt(np.pi)
    return result

def calibrate(cfg):
    r=cfg['risk']
    return brentq(lambda a:expected_risk(a,r['beta'],r['comorbidity_prevalence'])-r['prevalence'],-12,0,xtol=1e-13)

def fit_models(cfg):
    r=cfg['risk']; alpha=calibrate(cfg)
    rg=rng_for(cfg['seed'],'development')
    pop=patients(rg,(r['development_N'],),alpha,cfg)
    c=record_binary(pop.c,r['capture_A'],r['specificity'],rg.random(pop.c.shape))
    design=np.column_stack([np.ones(r['development_N']),pop.x,c])
    y=pop.outcome
    def loss(coef):
        eta=design@coef
        return np.mean(np.logaddexp(0,eta)-y*eta),design.T@(expit(eta)-y)/len(y)
    fit=minimize(loss,np.r_[alpha,r['beta']],jac=True,method='BFGS',options={'gtol':1e-9})
    if np.linalg.norm(fit.jac,np.inf)>1e-7:
        raise RuntimeError(f'Risk fit did not converge: {fit.message}')
    result={'alpha':alpha,'B1':list(np.r_[alpha,r['beta']]),'B2':fit.x.tolist(),
            'development_N':r['development_N'],'gradient_inf':float(np.max(np.abs(fit.jac))),
            'development_event_rate':float(y.mean()),'development_capture':r['capture_A'],
            'fit_message':str(fit.message)}
    write_json(result,'results/model_parameters.json')
    return result

def recorded_expected(coef,capture,specificity,cfg):
    pi=np.asarray(cfg['risk']['comorbidity_prevalence'])
    return expected_risk(coef[0],coef[1:],pi*capture+(1-pi)*(1-specificity))
