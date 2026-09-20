"""One-command pipeline, executed from project root."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import argparse,json,hashlib,time,platform,datetime
from src.common import ROOT,config,write_json
from src.risk_adjustment import fit_models
from src.simulations import run_main,run_sensitivity,run_reliability
from src.spc import run_spc

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=['all','main','sensitivity','spc','build'],default='all');parser.add_argument('--smoke',action='store_true');args=parser.parse_args()
    cfg=config()
    if args.smoke:
        cfg['replicates']=100;cfg['risk']['development_N']=10000
        raise SystemExit('Smoke runs use tests/test_smoke.py to avoid overwriting full-run results. Run python -m pytest tests/test_smoke.py.')
    start=time.time()
    if args.stage=='all':
        write_json({'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'config_sha256':hashlib.sha256((ROOT/'config/simulation_config.yaml').read_bytes()).hexdigest(),'protocol_sha256':hashlib.sha256((ROOT/'manuscript/protocol_for_preregistration.md').read_bytes()).hexdigest(),'seed':cfg['seed'],'R':cfg['replicates'],'python':platform.python_version()},'audit/run_manifest.json')
    if args.stage in ['all','main','sensitivity']:
        model_path=ROOT/'results/model_parameters.json'
        models=fit_models(cfg) if args.stage in ['all','main'] or not model_path.exists() else json.loads(model_path.read_text())
    if args.stage in ['all','main']:run_main(cfg,models)
    if args.stage in ['all','sensitivity']:run_sensitivity(cfg,models);run_reliability(cfg)
    if args.stage in ['all','spc']:run_spc(cfg)
    if args.stage in ['all','build']:
        from src.reporting import build
        build(cfg)
    print(f'Completed {args.stage} in {time.time()-start:.1f} seconds',flush=True)
    if args.stage=='all':
        manifest=json.loads((ROOT/'audit/run_manifest.json').read_text());manifest.update(completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),elapsed_seconds=time.time()-start)
        write_json(manifest,'audit/run_manifest.json')

if __name__=='__main__':main()
