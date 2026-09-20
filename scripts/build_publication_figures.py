"""Rebuild revised figures 1 and 3 from existing results; no simulations."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import pandas as pd
from src.common import ROOT,config
from src.reporting import build_framework,build_robustness
if __name__=='__main__':
    build_framework()
    build_robustness(pd.read_csv(ROOT/'results/main/ac.csv'),pd.read_csv(ROOT/'results/main/risk_paired.csv'),pd.read_csv(ROOT/'tables/table2_tipping_points.csv'))
