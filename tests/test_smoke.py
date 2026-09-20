from src.common import config
from src.risk_adjustment import calibrate
from src.simulations import risk_curves
import numpy as np

def test_risk_small_run_and_identity():
    cfg=config();cfg['replicates']=100
    a=calibrate(cfg);models={'alpha':a,'B1':[a,*cfg['risk']['beta']],'B2':[a,*cfg['risk']['beta']]}
    rows=risk_curves(cfg,models,500,captures=[.5,.7,.9])
    assert len(rows)==18
    for row in rows:
        if row['c_B']==.7 and row['g']==0:assert row['delta']==0
        if row['c_B']==.9 and row['g']==0:assert row['delta']<0
