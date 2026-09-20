import numpy as np

def safety_expectation(p,se,sp=1):
    return p*se+(1-p)*(1-sp)

def process_expectation(e,q,d,f,n):
    """Exact E(N/D | D>0); unconditional ratio undefined if D=0."""
    return e*q*d*n/(e*d+(1-e)*f)

def summarize(a,b,true_change,cfg,direction=-1,se=None):
    a,b = np.asarray(a),np.asarray(b)
    delta = b-a
    rc = np.divide(delta,a,out=np.full_like(delta,np.nan),where=a>0)
    valid = np.isfinite(rc)
    out={'R':len(a),'valid_R':int(valid.sum()),'undefined_R':int((~valid).sum()),
         'mean_A':float(np.nanmean(a)), 'mean_B':float(np.nanmean(b)),
         'delta':float(np.nanmean(delta)), 'bias':float(np.nanmean(delta)-true_change),
         'true_change':float(true_change),'delta_sd':float(np.nanstd(delta,ddof=1)),
         'delta_mcse':float(np.nanstd(delta,ddof=1)/np.sqrt(np.isfinite(delta).sum())),
         'mean_RC':float(np.nanmean(rc)),
         'ratio_of_means_RC':float(np.nanmean(b)/np.nanmean(a)-1),
         'RC_q025':float(np.nanquantile(rc,.025)),'RC_q975':float(np.nanquantile(rc,.975)),
         'delta_q025':float(np.nanquantile(delta,.025)), 'delta_q975':float(np.nanquantile(delta,.975))}
    def probability(name,event,denom=valid):
        k=int(np.sum(event & denom)); n=int(np.sum(denom)); p=k/n if n else np.nan
        out[name]=p; out[name+'_mcse']=float(np.sqrt(p*(1-p)/n)) if n else np.nan
        # Wilson bounds remain informative even when p=0 or 1.
        z=1.959963984540054
        mid=(p+z*z/(2*n))/(1+z*z/n)
        half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/(1+z*z/n)
        out[name+'_lo']=mid-half; out[name+'_hi']=mid+half
    for t in cfg['thresholds']:
        tag=str(round(t*100))
        probability('improve_'+tag,direction*rc>t)
        probability('deteriorate_'+tag,direction*rc < -t)
    probability('direction_reversal',direction*rc<0)
    probability('masked',np.abs(rc)<=cfg['mask_band'])
    if se is not None:
        out['mean_naive_CI_width']=float(np.nanmean(2*1.959963984540054*se))
        probability('naive_CI_excludes_zero',np.abs(delta)>1.959963984540054*se,np.isfinite(delta))
    return out
