import numpy as np
from .common import rng_for,save_csv

def evaluate_chart(rates,baseline,run_length=8,sigma=3,n=1000):
    centre=rates[:,:baseline].mean(axis=1)
    sd=np.sqrt(centre*(1-centre)/n)
    lower=np.maximum(0,centre-sigma*sd);upper=np.minimum(1,centre+sigma*sd)
    r=len(rates);first_run=np.full(r,np.nan);first_sigma=np.full(r,np.nan)
    pos=np.zeros(r,dtype=int);neg=pos.copy()
    for j in range(baseline,rates.shape[1]):
        value=rates[:,j]
        pos=np.where(value>centre,pos+1,0);neg=np.where(value<centre,neg+1,0)
        hit=(pos>=run_length)|(neg>=run_length)
        first_run=np.where(np.isnan(first_run)&hit,j+1,first_run)
        hit=(value<lower)|(value>upper)
        first_sigma=np.where(np.isnan(first_sigma)&hit,j+1,first_sigma)
    return centre,lower,upper,first_run,first_sigma

def run_spc(cfg):
    c=cfg['spc'];s=cfg['safety'];rows=[]
    for n in [c['monthly_N'],*c['monthly_N_sensitivity']]:
        for kind in ['abrupt','gradual']:
            for end in c['endpoints']:
                R=cfg['replicates'];rg=rng_for(cfg['seed'],'SPC',n,kind,end)
                sensitivity=np.r_[np.repeat(s['se_A'],c['baseline_months']),np.repeat(end,c['followup_months']) if kind=='abrupt' else np.linspace(s['se_A'],end,c['followup_months'])]
                truth=rg.binomial(n,s['p'],size=(R,len(sensitivity)))
                counts=rg.binomial(truth,sensitivity)
                rates=counts/n
                centre,low,high,fr,fs=evaluate_chart(rates,c['baseline_months'],c['run_length'],c['sigma'],n)
                for rule,first in [('eight_same_side',fr),('three_sigma',fs)]:
                    signal=np.isfinite(first);p=signal.mean();lag=first-(c['baseline_months']+1)
                    z=1.959963984540054;mid=(p+z*z/(2*R))/(1+z*z/R);half=z*np.sqrt(p*(1-p)/R+z*z/(4*R*R))/(1+z*z/R)
                    shift=rates[:,c['baseline_months']:].mean(axis=1)-centre
                    rows.append(dict(n=n,kind=kind,end=end,rule=rule,R=R,signal_probability=p,mcse=np.sqrt(p*(1-p)/R),wilson_lo=mid-half,wilson_hi=mid+half,detected_R=signal.sum(),nondetected_R=(~signal).sum(),mean_first_month_conditional=np.nanmean(first),median_lag_conditional=np.nanmedian(lag),mean_lag_conditional=np.nanmean(lag),restricted_mean_lag=np.mean(np.where(signal,lag,c['followup_months'])),mean_shift=shift.mean(),shift_mcse=shift.std(ddof=1)/np.sqrt(R)))
                if n==c['monthly_N'] and kind=='abrupt' and end==.9:
                    save_csv([dict(month=j+1,true_probability=s['p'],true_realised=truth[0,j]/n,sensitivity=sensitivity[j],observed=rates[0,j],centre=centre[0],lower=low[0],upper=high[0],first_run=fr[0],first_sigma=fs[0]) for j in range(len(sensitivity))],'results/main/spc_example.csv')
    save_csv(rows,'results/main/spc.csv')
