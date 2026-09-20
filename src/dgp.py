from dataclasses import dataclass
import numpy as np
from scipy.special import expit

@dataclass(frozen=True)
class Patients:
    x: np.ndarray
    c: np.ndarray
    risk: np.ndarray
    outcome: np.ndarray
    eligibility: np.ndarray
    action: np.ndarray

def patients(rng, shape, alpha, cfg):
    x = rng.normal(size=(*shape, 2))
    c = rng.random((*shape, 3)) < cfg['risk']['comorbidity_prevalence']
    beta = np.asarray(cfg['risk']['beta'])
    p = expit(alpha + x @ beta[:2] + c @ beta[2:])
    y = rng.random(shape) < p
    e = rng.random(shape) < cfg['process']['eligibility']
    a = e & (rng.random(shape) < cfg['process']['completion'])
    values = [x, c, p, y, e, a]
    for value in values:
        value.flags.writeable = False
    return Patients(*values)
