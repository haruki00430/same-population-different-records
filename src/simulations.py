import numpy as np
from scipy.special import expit
from .common import rng_for, grid, save_csv
from .recording import paired_counts, process_joint, record_binary
from .measures import summarize, safety_expectation, process_expectation
from .risk_adjustment import recorded_expected

def safety_cell(cfg,N,p,se_a,se_b,g=0,sp_a=1,sp_b=1,design='paired',coupling='shared',family='main'):
    R=cfg['replicates']; rg=rng_for(cfg['seed'],'A',N,p,se_a,g,sp_a,sp_b,design,coupling,family)
    pb=p*(1-g)
    if design=='paired':
        # Latent states: events at both times, events prevented, non-events.
        latent=rg.multinomial(N,[pb,p-pb,1-p],size=R)
        ca=np.zeros(R,dtype=int); cb=ca.copy(); discord=np.zeros(R,dtype=int)
        for j,aa,bb in [(0,se_a,se_b),(1,se_a,1-sp_b),(2,1-sp_a,1-sp_b)]:
            a,b=paired_counts(rg,latent[:,j],aa,bb,coupling=='independent')
            ca+=a; cb+=b
            # Discordance variance calculated exactly from joint draw when coupling shared;
            # independent recording needs joint counts, regenerated below via same helper logic.
            # For shared uniform, nested events guarantee |a-b| within each latent stratum.
            if coupling=='shared': discord+=np.abs(a-b)
        if coupling=='independent':
            # Re-run with explicit 4-category counts to retain discordance.
            rg=rng_for(cfg['seed'],'A-independent',N,p,se_a,se_b,g,sp_a,sp_b,family)
            ca[:]=0;cb[:]=0;discord[:]=0
            for j,aa,bb in [(0,se_a,se_b),(1,se_a,1-sp_b),(2,1-sp_a,1-sp_b)]:
                d=rg.multinomial(latent[:,j],[aa*bb,aa*(1-bb),(1-aa)*bb,(1-aa)*(1-bb)])
                ca+=d[:,0]+d[:,1];cb+=d[:,0]+d[:,2];discord+=d[:,1]+d[:,2]
        a=ca/N;b=cb/N
        se=np.sqrt(np.maximum(discord/N-(b-a)**2,0)/(N-1))
    else:
        a=rg.binomial(N,safety_expectation(p,se_a,sp_a),R)/N
        b=rg.binomial(N,safety_expectation(pb,se_b,sp_b),R)/N
        se=np.sqrt(a*(1-a)/N+b*(1-b)/N)
    out=dict(architecture='A',family=family,design=design,coupling=coupling,N=N,p=p,g=g,se_A=se_a,se_B=se_b,sp_A=sp_a,sp_B=sp_b)
    out.update(summarize(a,b,-p*g,cfg,se=se))
    ea=safety_expectation(p,se_a,sp_a);eb=safety_expectation(pb,se_b,sp_b)
    out.update(expected_A=ea,expected_B=eb,expected_delta=eb-ea,expected_RC=eb/ea-1)
    out['analytic_z_delta']=(out['delta']-(eb-ea))/out['delta_mcse'] if out['delta_mcse']>0 else 0.
    return out

def process_cell(cfg,N,n_b=None,d_b=None,f_b=None,g=0,e=None,q=None,n_a=None,design='paired',coupling='shared',family='C1'):
    c=cfg['process'];e=c['eligibility'] if e is None else e;q=c['completion'] if q is None else q
    na=c['se_N'] if n_a is None else n_a;nb=na if n_b is None else n_b
    da=c['se_D'];db=da if d_b is None else d_b;fa=c['fpr_D'];fb=fa if f_b is None else f_b
    qb=q*(1+g)
    if qb>1: raise ValueError('Completion exceeds 1')
    rg=rng_for(cfg['seed'],'C',N,e,q,na,g,design,coupling,family)
    joint=process_joint(e,q,qb,da,db,fa,fb,na,nb,coupling=='independent')
    R=cfg['replicates']
    if design=='paired':
        counts=rg.multinomial(N,joint.ravel(),size=R).reshape(R,3,3)
        aa=counts.sum(axis=2);bb=counts.sum(axis=1)
    else:
        aa=rg.multinomial(N,joint.sum(axis=1),size=R)
        bb=rg.multinomial(N,joint.sum(axis=0),size=R)
    dena=aa[:,1:].sum(axis=1);denb=bb[:,1:].sum(axis=1)
    a=np.divide(aa[:,2],dena,out=np.full(R,np.nan),where=dena>0)
    b=np.divide(bb[:,2],denb,out=np.full(R,np.nan),where=denb>0)
    out=dict(architecture='C',family=family,design=design,coupling=coupling,N=N,e=e,q=q,g=g,n_A=na,n_B=nb,d_A=da,d_B=db,f_A=fa,f_B=fb)
    out.update(summarize(a,b,q*g,cfg,direction=1))
    ea=process_expectation(e,q,da,fa,na);eb=process_expectation(e,qb,db,fb,nb)
    out.update(expected_A=ea,expected_B=eb,expected_delta=eb-ea,expected_RC=eb/ea-1)
    out['analytic_z_delta']=(out['delta']-(eb-ea))/out['delta_mcse'] if out['delta_mcse']>0 else 0.
    return out

def risk_curves(cfg,models,N,design='paired',family='main',capture_a=None,captures=None,mode='independent',specificity=None):
    r=cfg['risk']; R=cfg['replicates'];batch=cfg['batch_size'];ca=r['capture_A'] if capture_a is None else capture_a
    captures=grid(r['capture_grid']) if captures is None else np.asarray(captures)
    sp=r['specificity'] if specificity is None else specificity
    gs=cfg['improvements'] if family=='main' and design=='paired' else [0.]
    a_store=np.zeros((2,len(gs),R));b_store=np.zeros((2,len(gs),len(captures),R))
    rg=rng_for(cfg['seed'],'B',N,design,family,ca,mode,sp)
    beta=np.array(r['beta']);coefs=[np.array(models[k]) for k in ['B1','B2']]
    def sample(m):
        x=rg.normal(size=(m,N,2));c=rg.random((m,N,3))<r['comorbidity_prevalence']
        p=expit(models['alpha']+x@beta[:2]+c@beta[2:]);u=rg.random((m,N))
        uc=rg.random((m,N,3))
        if mode=='correlated':
            uc=np.where(rg.random((m,N,1))<r['correlated_error_mixture'],rg.random((m,N,1)),uc)
        return x,c,p,u,uc
    def captured(c,u,level):
        sens=level
        if mode=='severity':sens=np.clip(level+np.where(c[...,2:3],r['severity_capture_offset'],-r['severity_capture_offset']),0,1)
        return record_binary(c,sens,sp,u)
    for start in range(0,R,batch):
        stop=min(start+batch,R);m=stop-start
        xa,c,p,u,uc=sample(m)
        coded=captured(c,uc,ca)
        xb,cb,pb,ub,ucb=(xa,c,p,u,uc) if design=='paired' else sample(m)
        if mode=='independent_records': ucb=rg.random(uc.shape)
        obs_a=(u<p).sum(axis=1)
        obs_b=[(ub<pb*(1-g)).sum(axis=1) for g in gs]
        base_linear=[co[0]+xa@co[1:3] for co in coefs]
        b_linear=[co[0]+xb@co[1:3] for co in coefs]
        for k,co in enumerate(coefs):
            expected=expit(base_linear[k]+coded@co[3:]).sum(axis=1)
            for h in range(len(gs)):a_store[k,h,start:stop]=obs_a/expected
        for j,level in enumerate(captures):
            coded_b=captured(cb,ucb,level)
            for k,co in enumerate(coefs):
                expected=expit(b_linear[k]+coded_b@co[3:]).sum(axis=1)
                for h in range(len(gs)):b_store[k,h,j,start:stop]=obs_b[h]/expected
    rows=[]
    for k,name in enumerate(['B1','B2']):
        ea=recorded_expected(coefs[k],ca,sp,cfg)
        for h,g in enumerate(gs):
            for j,level in enumerate(captures):
                eb=recorded_expected(coefs[k],level,sp,cfg)
                out=dict(architecture=name,family=family,design=design,N=N,g=g,c_A=ca,c_B=float(level),specificity=sp,mode=mode)
                out.update(summarize(a_store[k,h],b_store[k,h,j],-g,cfg))
                if mode not in ['correlated','severity']:
                    out.update(population_OE_A=r['prevalence']/ea,population_OE_B=r['prevalence']*(1-g)/eb,expected_RC=(1-g)*ea/eb-1)
                rows.append(out)
    return rows

def run_main(cfg,models):
    n=cfg['N'];s=cfg['safety'];c=cfg['process'];rows=[]
    for design in ['paired','repeated']:
        for g in (cfg['improvements'] if design=='paired' else [0.]):
            for se in grid(s['se_grid']):rows.append(safety_cell(cfg,n,s['p'],s['se_A'],float(se),g,design=design))
            for v in grid(c['numerator_grid']):rows.append(process_cell(cfg,n,n_b=float(v),g=g,design=design))
        for v in grid(c['denominator_grid']):rows.append(process_cell(cfg,n,d_b=float(v),design=design,family='C2_sensitivity'))
        for v in grid(c['fpr_grid']):rows.append(process_cell(cfg,n,f_b=float(v),design=design,family='C2_false_positive'))
    for v in grid(c['numerator_grid']):
        for d in grid(c['denominator_grid']):rows.append(process_cell(cfg,n,n_b=float(v),d_b=float(d),family='C3_surface'))
    save_csv(rows,'results/main/ac.csv')
    for design in ['paired','repeated']:
        print('Risk main',design,flush=True)
        save_csv(risk_curves(cfg,models,n,design),f'results/main/risk_{design}.csv')

def run_sensitivity(cfg,models):
    n=cfg['N'];s=cfg['safety'];c=cfg['process'];rows=[]
    variants=[('sample_size',{'N':v}) for v in cfg['sample_sizes']]
    variants += [('prevalence',{'p':v}) for v in s['prevalence_sensitivity']]
    variants += [('specificity_equal',{'sp_a':v,'sp_b':v}) for v in s['specificity_sensitivity']]
    variants += [('specificity_drift',{'sp_b':v}) for v in s['specificity_sensitivity']]
    variants += [('baseline_capture',{'se_a':v}) for v in s['baseline_sensitivity']]
    variants += [('independent_recording',{'coupling':'independent'})]
    for family,overrides in variants:
        base=dict(N=n,p=s['p'],se_a=s['se_A']);base.update(overrides)
        for design in ['paired','repeated']:
            for se in grid(s['se_grid']):rows.append(safety_cell(cfg,**base,se_b=float(se),design=design,family=family))
    variants=[('sample_size',{'N':v}) for v in cfg['sample_sizes']]
    variants += [('eligibility',{'e':v}) for v in c['eligibility_sensitivity']]
    variants += [('completion',{'q':v}) for v in c['completion_sensitivity']]
    variants += [('baseline_capture',{'n_a':v}) for v in [.6,.9]]
    variants += [('independent_recording',{'coupling':'independent'})]
    for family,overrides in variants:
        base=dict(N=n);base.update(overrides)
        for v in grid(c['numerator_grid']):rows.append(process_cell(cfg,**base,n_b=float(v),family=family))
    # Negative control: eliminate denominator false positives in BOTH records.
    import copy
    perfect=copy.deepcopy(cfg);perfect['process']['fpr_D']=0
    for d in grid(c['denominator_grid']):rows.append(process_cell(perfect,n,d_b=float(d),family='perfect_denominator_specificity'))
    save_csv(rows,'results/sensitivity/ac.csv')
    risk_rows=[];r=cfg['risk']
    variants=[('sample_size',{'N':v}) for v in cfg['sample_sizes']]
    variants += [('baseline_capture',{'capture_a':v}) for v in r['baseline_sensitivity']]
    variants += [('correlated',{'mode':'correlated'}),('severity',{'mode':'severity'}),('independent_recording',{'mode':'independent_records'})]
    variants += [('specificity',{'specificity':v}) for v in [1.,.99]]
    for family,overrides in variants:
        base=dict(N=n);base.update(overrides)
        print('Risk sensitivity',family,overrides,flush=True)
        risk_rows+=risk_curves(cfg,models,**base,family=family,captures=r['sensitivity_grid'])
        save_csv(risk_rows,'results/sensitivity/risk.csv')

def run_reliability(cfg):
    rows=[];lo,hi=cfg['reliability']['prevalence_range'];sea=cfg['safety']['se_A'];seb=cfg['reliability']['se_B']
    for N in [500,1000,2000,5000,10000]:
        rg=rng_for(cfg['seed'],'reliability',N);p=rg.uniform(lo,hi,cfg['replicates'])
        latent=rg.binomial(N,p);a,b=paired_counts(rg,latent,sea,seb)
        for label,se,obs in [('A',sea,a/N),('B',seb,b/N)]:
            between=se**2*(hi-lo)**2/12
            mean=se*(lo+hi)/2; second=se**2*(lo*lo+lo*hi+hi*hi)/3
            noise=(mean-second)/N
            est_noise=np.mean(obs*(1-obs)/(N-1));est_between=max(np.var(obs,ddof=1)-est_noise,0)
            rows.append(dict(N=N,record=label,R=cfg['replicates'],analytic_reliability=between/(between+noise),estimated_reliability=est_between/(est_between+est_noise),mean_measure=obs.mean(),expected_measure=mean,expected_RC=seb/sea-1,mean_paired_delta=np.mean((b-a)/N)))
    save_csv(rows,'results/sensitivity/reliability.csv')
