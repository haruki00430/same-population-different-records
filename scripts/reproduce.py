import subprocess,sys
from pathlib import Path
root=Path(__file__).resolve().parents[1]
def run(*args):subprocess.run([sys.executable,*args],cwd=root,check=True)
run('-m','pytest','-q','tests/test_core.py','tests/test_smoke.py')
run('scripts/run_all.py','--stage','all')
with (root/'audit/test_results.txt').open('w',encoding='utf-8') as out:
    subprocess.run([sys.executable,'-m','pytest','-q'],cwd=root,stdout=out,stderr=subprocess.STDOUT,check=True)
run('scripts/build_manuscript.py')
