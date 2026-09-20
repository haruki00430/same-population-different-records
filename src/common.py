from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]

def config():
    return yaml.safe_load((ROOT / 'config/simulation_config.yaml').read_text())

def rng_for(seed, *labels):
    # No dependence on Python's randomized hash or execution order.
    digest = hashlib.sha256(json.dumps(labels, sort_keys=True).encode()).digest()
    return np.random.default_rng(np.random.SeedSequence([seed, *np.frombuffer(digest[:16], dtype='<u4')]))

def grid(spec):
    lo, hi, step = spec
    return np.round(np.linspace(lo, hi, round((hi-lo)/step)+1), 8)

def save_csv(rows, path):
    path = ROOT / path
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False, float_format='%.12g')

def write_json(obj, path):
    p = ROOT / path
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(obj, indent=2), encoding='utf-8')
