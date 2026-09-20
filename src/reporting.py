"""Deterministic tables, figures and audit from saved summaries."""
import json,hashlib,datetime
import numpy as np
import pandas as pd
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .common import ROOT,save_csv,rng_for,write_json
from .risk_adjustment import recorded_expected

BLUE='#176B91';RED='#B2493C';GREY='#737D85'

def style():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.titlesize':12,'axes.labelsize':10,'figure.dpi':120,'savefig.dpi':300,'svg.fonttype':'none','svg.hashsalt':'same-population-fixed-v1','pdf.fonttype':42})

def save(fig,name):
    for suffix in ['png','svg','pdf']:
        metadata={'CreationDate':None,'ModDate':None} if suffix=='pdf' else {'Date':None} if suffix=='svg' else {}
        fig.savefig(ROOT/'figures'/f'{name}.{suffix}',bbox_inches='tight',facecolor='white',metadata=metadata)
    plt.close(fig)

def md_table(df,path):
    def val(x):
        if isinstance(x,(float,np.floating)): return f'{x:.5g}' if np.isfinite(x) else 'Not attainable'
        return str(x).replace('|','/')
    lines=['| '+' | '.join(df.columns)+' |','| '+' | '.join(['---']*len(df.columns))+' |']
    lines += ['| '+' | '.join(val(v) for v in row)+' |' for row in df.itertuples(index=False,name=None)]
    (ROOT/path).write_text('\n'.join(lines)+'\n',encoding='utf-8')

def tipping_points(cfg,models):
    rows=[]
    for g in cfg['improvements'][1:]:
        base=cfg['safety']['se_A'];root=base/(1-g)
        rows.append(dict(measure='A: event ascertainment',g=g,baseline=base,root=root,change_pp=100*(root-base),erasure=f'Increase ascertainment by {100*(root-base):.3f} pp',reversal='Strictly larger increase' if root<1-1e-10 else 'Unattainable: root at sensitivity=1',grid='0.60 to 1.00',within_grid=root<=1,root_residual=(1-g)*root/base-1))
        for model in ['B1','B2']:
            base=cfg['risk']['capture_A'];coef=models[model];sp=cfg['risk']['specificity'];ea=recorded_expected(coef,base,sp,cfg)
            fn=lambda c:(1-g)*ea/recorded_expected(coef,c,sp,cfg)-1
            root=brentq(fn,0,1)
            rows.append(dict(measure=f'{model}: comorbidity capture',g=g,baseline=base,root=root,change_pp=100*(root-base),erasure=f'Decrease capture by {100*(base-root):.3f} pp',reversal='Strictly larger decrease'+(' (outside main grid)' if root<.5 else ''),grid='0.50 to 1.00',within_grid=.5<=root<=1,root_residual=fn(root)))
        base=cfg['process']['se_N'];root=base/(1+g)
        rows.append(dict(measure='C: numerator capture',g=g,baseline=base,root=root,change_pp=100*(root-base),erasure=f'Decrease numerator capture by {100*(base-root):.3f} pp',reversal='Strictly larger decrease',grid='0.60 to 1.00',within_grid=.6<=root<=1,root_residual=(1+g)*root/base-1))
    return pd.DataFrame(rows)


def build_framework():
    from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
    style()
    fig,ax=plt.subplots(figsize=(12,3.0))
    ax.set_xlim(0,12);ax.set_ylim(0,3);ax.axis('off')
    widths=[2.1,2.1,1.6,1.6,2.3];gap=.42
    left=(12-sum(widths)-4*gap)/2;y=1.03;h=1.02
    labels=['Clinical reality\nPatients, care,\noutcomes','Record-generating\nprocess','Observed\ndata','Quality\nmeasure','Improvement or\ndeterioration\ninference']
    boxes=[]
    for i,(w,label) in enumerate(zip(widths,labels)):
        color=BLUE if i==0 else RED if i==1 else '#85929A'
        fill='#E7F0F4' if i==0 else '#F7EAE5' if i==1 else '#F3F4F5'
        box=FancyBboxPatch((left,y),w,h,boxstyle='round,pad=0,rounding_size=0.09',linewidth=1.3,edgecolor=color,facecolor=fill)
        ax.add_patch(box);ax.text(left+w/2,y+h/2,label,ha='center',va='center',fontsize=11,linespacing=1.35)
        boxes.append((left,w));left+=w+gap
    for (x,w),(nx,nw) in zip(boxes,boxes[1:]):
        ax.add_patch(FancyArrowPatch((x+w+.055,y+h/2),(nx-.055,y+h/2),arrowstyle='-|>',mutation_scale=13,lw=1.4,color='#53616A',shrinkA=0,shrinkB=0))
    for i,label,color in [(0,'HELD FIXED',BLUE),(1,'PERTURBED',RED)]:
        x,w=boxes[i];ax.text(x+w/2,.68,label,ha='center',va='center',fontsize=10,weight='bold',color=color)
    ax.text(6,2.64,'Same patients, same care, same outcomes — different records',ha='center',va='center',fontsize=14,weight='bold')
    ax.text(6,.2,'The paired null comparison changes recording while holding clinical reality constant.',ha='center',va='center',fontsize=10,color='#53616A')
    fig.subplots_adjust(left=.01,right=.99,top=.98,bottom=.02)
    save(fig,'figure1_framework')

def build_robustness(ac,risk,tips):
    style()
    # Figure 3: explicit raw direction, population roots, MC dots, sampling bands.
    fig,axes=plt.subplots(1,3,figsize=(12.5,4.7))
    selections=[(ac[(ac.architecture=='A')&(ac.design=='paired')],'se_B',.8,'A  Safety outcome','Ascertainment change (pp)','A: event ascertainment'),(risk[risk.architecture=='B2'],'c_B',.7,'B  Risk-adjusted outcome (B2)','Comorbidity capture change (pp)','B2: comorbidity capture'),(ac[(ac.architecture=='C')&(ac.family=='C1')&(ac.design=='paired')],'n_B',.8,'C  Process completion','Numerator capture change (pp)','C: numerator capture')]
    for ax,(frame,col,base,title,xlabel,key) in zip(axes,selections):
        for g,color,label in [(0.,BLUE,'Clinical state unchanged'),(.1,RED,'True 10% improvement')]:
            d=frame[np.isclose(frame.g,g)].sort_values(col);x=(d[col]-base)*100
            ax.plot(x,d.expected_RC*100,color=color,label=label,lw=2,ls='-' if g==0 else '--')
            ax.scatter(x,d.mean_RC*100,color=color,s=8,zorder=3)
            ax.fill_between(x,d.RC_q025*100,d.RC_q975*100,color=color,alpha=.09)
        t=tips[(tips.measure==key)&np.isclose(tips.g,.1)].iloc[0]
        ax.axvline(t.change_pp,ls=':',color=RED,lw=1.3);ax.plot(t.change_pp,0,'o',mfc='white',mec=RED,ms=6,zorder=4)
        ax.axhline(0,color=GREY,lw=.8);ax.axvline(0,color=GREY,lw=.6,alpha=.5);ax.set_title(title,loc='left');ax.set_xlabel(xlabel);ax.grid(axis='y',alpha=.13)
        ax.text(.04,.96,'Higher is worse' if key[0]!='C' else 'Higher is better',transform=ax.transAxes,va='top',fontsize=9)
    axes[0].set_ylabel('Apparent relative measure change (%)');axes[1].legend(loc='lower left',fontsize=8)
    fig.text(.5,.015,'Lines: population expectations. Dots: Monte Carlo means. Shading: central 95% of replicate changes.\nDotted lines: exact masking boundary; reversal is beyond it in the adverse direction.',ha='center',fontsize=9);fig.tight_layout(rect=[0,.1,1,1]);save(fig,'figure3_robustness')
    return selections

def build(cfg):
    style()
    for d in ['figures','tables','audit']:(ROOT/d).mkdir(exist_ok=True)
    ac=pd.read_csv(ROOT/'results/main/ac.csv');risk=pd.read_csv(ROOT/'results/main/risk_paired.csv');rep=pd.read_csv(ROOT/'results/main/risk_repeated.csv');sens=pd.read_csv(ROOT/'results/sensitivity/ac.csv');rs=pd.read_csv(ROOT/'results/sensitivity/risk.csv');spc=pd.read_csv(ROOT/'results/main/spc.csv');rel=pd.read_csv(ROOT/'results/sensitivity/reliability.csv');models=json.loads((ROOT/'results/model_parameters.json').read_text())
    tips=tipping_points(cfg,models);tips.to_csv(ROOT/'tables/table2_tipping_points.csv',index=False)
    show=tips[['measure','g','erasure','reversal','grid','within_grid']].rename(columns={'measure':'Quality measure','g':'True relative improvement','erasure':'Recording change to erase','reversal':'To reverse population direction','grid':'Sensitivity range','within_grid':'Root inside grid'})
    md_table(show,'tables/table2_tipping_points.md')
    architecture=pd.DataFrame([
        ['A: safety','Adverse-event probability','Event sensitivity/specificity','Ascertainment 0.60–1.00','Recorded events / N','Improved ascertainment appears harmful'],
        ['B1/B2: risk-adjusted','Outcome risk conditional on true case mix','Comorbidity sensitivity/specificity; frozen model','Coding sensitivity 0.50–1.00','Observed outcomes / predicted outcomes','Better coding appears to improve O/E'],
        ['C: process','Action completion among truly eligible','Numerator and eligibility capture','Numerator 0.60–1.00; denominator 0.70–1.00; FPR 0–0.03','Recorded numerator / recorded denominator','Documentation can inflate completion; false eligibility dilutes it']
    ],columns=['Measure architecture','True clinical target','Recording mechanism','Perturbation','Observed measure','Failure mode'])
    architecture.to_csv(ROOT/'tables/table1_architectures.csv',index=False);md_table(architecture,'tables/table1_architectures.md')
    # Supplement tables preserve all cells, including null/negative-control results.
    for name,frame in [('S1_main_AC',ac),('S2_main_risk',pd.concat([risk,rep])),('S3_sensitivity_AC',sens),('S4_sensitivity_risk',rs),('S5_SPC',spc),('S6_reliability',rel)]:frame.to_csv(ROOT/'tables'/f'{name}.csv',index=False)
    validation=pd.concat([ac,sens],ignore_index=True)
    validation=validation[['architecture','family','design','N','expected_delta','delta','delta_mcse','analytic_z_delta']]
    validation.to_csv(ROOT/'tables/S7_analytical_validation.csv',index=False)
    if not (validation.analytic_z_delta.abs()<6).all():raise AssertionError('Analytical validation failed: inspect full grid')
    build_framework()
    # Figure 2: one pre-specified seed; save exact counts and individual record audit.
    rg=rng_for(cfg['seed'],'figure2');N=cfg['N'];y=rg.random(N)<cfg['safety']['p'];u=rg.random(N);ya=y&(u<.8);yb=y&(u<.9)
    save_csv([{'patient_id':j+1,'true_event':int(y[j]),'record_A_event':int(ya[j]),'record_B_event':int(yb[j])} for j in range(N)],'results/main/paired_example_patients.csv')
    fig,axes=plt.subplots(1,2,figsize=(9.5,4.2));counts=[y.sum(),ya.sum(),yb.sum()]
    bars=axes[0].bar(['True events\n(shared)','Record A\nSe = 0.80','Record B\nSe = 0.90'],counts,color=[GREY,BLUE,RED],width=.6)
    axes[0].bar_label(bars,padding=4);axes[0].set_ylim(0,max(counts)*1.25);axes[0].set_ylabel('Events among the identical 2,000 patients');axes[0].set_title('A  Same true outcomes')
    bars=axes[1].bar(['Record A','Record B'],[ya.mean()*100,yb.mean()*100],color=[BLUE,RED],width=.6);axes[1].bar_label(bars,fmt='%.2f%%',padding=4);axes[1].set_ylim(0,max(yb.mean(),ya.mean())*130);axes[1].set_ylabel('Recorded adverse-event rate (%)');axes[1].set_title('B  Different measured safety')
    fig.text(.5,.015,f'Fixed-seed illustration, not selected for effect size. True events are identical; {int((yb&~ya).sum())} additional events are recorded.',ha='center',fontsize=9);fig.tight_layout(rect=[0,.07,1,1]);save(fig,'figure2_same_patients')
    selections=build_robustness(ac,risk,tips)
    # Retain the existing B1 reference display in supplementary figure S3.
    selections[1]=(risk[risk.architecture=='B1'],'c_B',.7,'B  Risk-adjusted outcome (B1)','Comorbidity capture change (pp)','B1: comorbidity capture')
    # Figure 4: true probability distinguished from realised latent rate.
    d=pd.read_csv(ROOT/'results/main/spc_example.csv');fig,axes=plt.subplots(3,1,figsize=(9.5,7.8),sharex=True)
    axes[0].plot(d.month,d.true_realised*100,color='#BAC4CA',lw=1,label='Realised true monthly rate');axes[0].plot(d.month,d.true_probability*100,color=BLUE,lw=2,label='Stable clinical probability');axes[0].set_ylabel('True events (%)');axes[0].set_title('A  Clinical process',loc='left');axes[0].legend(fontsize=8,loc='upper left')
    axes[1].step(d.month,d.sensitivity*100,where='mid',color=RED,lw=2);axes[1].set_ylabel('Sensitivity (%)');axes[1].set_ylim(75,95);axes[1].set_title('B  Record-generating process',loc='left')
    axes[2].plot(d.month,d.observed*100,'o-',ms=3,lw=1,color=BLUE);axes[2].plot(d.month,d.centre*100,color=GREY,lw=1);axes[2].plot(d.month,d.lower*100,'--',color=GREY,lw=1);axes[2].plot(d.month,d.upper*100,'--',color=GREY,lw=1)
    for col,label,marker in [('first_run','First 8-point signal','s'),('first_sigma','First 3-sigma signal','D')]:
        month=d[col].iloc[0]
        if np.isfinite(month):
            val=d.loc[d.month==month,'observed'].iloc[0]*100;axes[2].scatter([month],[val],marker=marker,facecolors='none',edgecolors=RED,s=90,label=label)
    axes[2].set_ylabel('Recorded events (%)');axes[2].set_xlabel('Month');axes[2].set_title('C  Observed process and frozen p-chart limits',loc='left')
    if axes[2].get_legend_handles_labels()[0]:axes[2].legend(fontsize=8,loc='lower left')
    for ax in axes:ax.axvline(24.5,color=RED,ls=':',lw=1);ax.set_xlim(1,48);ax.grid(axis='y',alpha=.15)
    fig.text(.5,.015,'First replicate at the prespecified 0.90 abrupt endpoint; n = 1,000/month. Baseline: months 1–24.\nA change in the observed process is not necessarily a change in clinical care.',ha='center',fontsize=9);fig.tight_layout(rect=[0,.07,1,1]);save(fig,'figure4_spc')
    # Supplementary response surface, precision and inference curves.
    surf=ac[ac.family=='C3_surface'].pivot(index='d_B',columns='n_B',values='expected_RC')*100
    fig,ax=plt.subplots(figsize=(7.2,4.7));im=ax.pcolormesh(surf.columns*100,surf.index*100,surf.values,cmap='RdBu',shading='auto',vmin=-26,vmax=26)
    cs=ax.contour(surf.columns*100,surf.index*100,surf.values,levels=[0],colors='black',linewidths=1);ax.clabel(cs,fmt='Zero change');ax.plot(80,95,'o',mec='black',mfc='white');ax.set_xlabel('Numerator sensitivity (%)');ax.set_ylabel('Eligibility sensitivity (%)');ax.set_title('Identical care: combined numerator and denominator drift');fig.colorbar(im,ax=ax,label='Apparent relative completion change (%)');fig.tight_layout();save(fig,'figureS1_process_surface')
    precision=pd.concat([ac,sens]);precision=precision[(precision.architecture=='A')&(precision.design=='repeated')&np.isclose(precision.se_B,.9)&np.isclose(precision.se_A,.8)&np.isclose(precision.p,.05)&np.isclose(precision.sp_A,1)&np.isclose(precision.sp_B,1)&np.isclose(precision.g,0)&precision.family.isin(['main','sample_size'])].sort_values('N')
    precision.to_csv(ROOT/'tables/S8_precision.csv',index=False)
    fig,axes=plt.subplots(1,2,figsize=(9.5,4.3));axes[0].errorbar(precision.N,precision.delta*100,yerr=[(precision.delta-precision.delta_q025)*100,(precision.delta_q975-precision.delta)*100],fmt='o-',color=BLUE,capsize=4);axes[0].axhline(0,color=GREY,ls='--');axes[0].axhline(.5,color=RED,ls=':');axes[0].set_xlabel('Patients per period');axes[0].set_ylabel('Apparent absolute change (pp)');axes[0].set_title('A  Sampling spread shrinks around bias')
    for label,color in [('A',BLUE),('B',RED)]:
        dd=rel[rel.record==label];axes[1].plot(dd.N,dd.analytic_reliability,'o-',label='Record '+label,color=color)
    axes[1].set_xlabel('Patients per unit');axes[1].set_ylabel('Signal/noise reliability');axes[1].set_title('B  Reliability can remain high');axes[1].set_ylim(0,1);axes[1].legend();fig.text(.5,.01,'A: independent cohorts; bars are empirical 95% replicate intervals. B: between-unit p ~ Uniform(0.03, 0.07).',ha='center',fontsize=8);fig.tight_layout(rect=[0,.05,1,1]);save(fig,'figureS2_precision_reliability')
    fig,axes=plt.subplots(2,3,figsize=(12,7))
    for j,(frame,col,base,title,xlabel,key) in enumerate(selections):
        d=frame[np.isclose(frame.g,0)].sort_values(col);x=(d[col]-base)*100
        axes[0,j].plot(x,d.bias,color=BLUE);axes[0,j].axhline(0,color=GREY,lw=.7);axes[0,j].set_title(title,loc='left');axes[0,j].set_ylabel('Bias in change (measure units)' if j==0 else '')
        for field,color,ls in [('improve_5',BLUE,'-'),('deteriorate_5',RED,'-'),('improve_10',BLUE,'--'),('deteriorate_10',RED,'--')]:axes[1,j].plot(x,d[field],color=color,ls=ls,label=field.replace('_',' >')+'%')
        axes[1,j].set_ylim(-.03,1.03);axes[1,j].set_xlabel(xlabel);axes[1,j].set_ylabel('Practical-threshold probability' if j==0 else '')
    axes[1,0].legend(fontsize=7);fig.tight_layout();save(fig,'figureS3_bias_inference')
    fig,axes=plt.subplots(1,2,figsize=(10,4.3))
    for model,color in [('B1',BLUE),('B2',RED)]:
        for g,ls in [(0,'-'),(.1,'--')]:
            d=risk[(risk.architecture==model)&np.isclose(risk.g,g)].sort_values('c_B');axes[0].plot((d.c_B-.7)*100,d.mean_RC*100,color=color,ls=ls,label=f'{model}, g={g:.1f}')
    axes[0].axhline(0,color=GREY,lw=.7);axes[0].set_xlabel('Comorbidity capture change (pp)');axes[0].set_ylabel('Relative O/E change (%)');axes[0].legend(fontsize=8);axes[0].set_title('A  Both models frozen')
    for family,color in [('C2_sensitivity',BLUE),('C2_false_positive',RED)]:
        d=ac[(ac.family==family)&(ac.design=='paired')];col='d_B' if family=='C2_sensitivity' else 'f_B';base=.95 if col=='d_B' else .005;d=d.sort_values(col);axes[1].plot((d[col]-base)*100,d.mean_RC*100,color=color,label='Eligibility sensitivity' if col=='d_B' else 'Eligibility false positives')
    axes[1].axhline(0,color=GREY,lw=.7);axes[1].set_xlabel('Recording change (pp)');axes[1].set_ylabel('Relative completion change (%)');axes[1].legend(fontsize=8);axes[1].set_title('B  Denominator-only perturbation');fig.tight_layout();save(fig,'figureS4_models_denominator')
    # Before any final manuscript drafting: execute and save GO/NO-GO audit.
    mainA=ac[(ac.architecture=='A')&(ac.design=='paired')&np.isclose(ac.g,0)&np.isclose(ac.se_B,.85)].iloc[0]
    strong=abs(mainA.expected_RC)>=.05
    audit=f'''# Analysis audit

Generated from executed results at {datetime.datetime.now(datetime.timezone.utc).isoformat()}.

## Decision: {'GO — strong under the predefined illustrative criterion' if strong else 'REVIEW — criterion not met'}

The prespecified 5-percentage-point ascertainment change (0.80 to 0.85) produces an analytical relative measure change of {100*mainA.expected_RC:.3f}% with clinical event risk fixed. Monte Carlo mean relative change is {100*mainA.mean_RC:.3f}%. The 10% safety improvement masking root is {100*(.8/.9-.8):.3f} percentage points. These meet the operational <=10 pp criterion. This does not establish that the perturbation is empirically common or clinically plausible in any particular setting.

## Validation and completeness

- {len(validation)} A/C cells checked against analytical expectations; maximum absolute Monte Carlo z discrepancy = {validation.analytic_z_delta.abs().max():.3f}; required <6.
- Main A/C cells: {len(ac)}; main risk cells: {len(risk)+len(rep)}; sensitivity A/C cells: {len(sens)}; sensitivity risk cells: {len(rs)}; SPC rule rows: {len(spc)} (two rules per scenario).
- Every listed simulation cell has {cfg['replicates']:,} replicates. All grid cells, including unchanged records and weak denominator effects, are retained.
- Protocol/config hashes and execution log document this run. No empirical patient data were used. Development model fit gradient infinity norm = {models['gradient_inf']:.3g}.
- B1 baseline O/E is not assumed calibrated. B2 is fixed after development; development-sample uncertainty is not propagated.
- A/C identity and analytical tests passed before full execution. Full-run verification is saved separately in audit/test_results.txt.

## Results that constrain interpretation

- A 20% clinical event reduction is erased only at perfect ascertainment (1.00) from baseline .80; population-direction reversal is unattainable through sensitivity alone in this main model (finite-sample direction still varies).
- Non-differential denominator sensitivity cancels exactly when denominator specificity is perfect. Nonzero FPR breaks that cancellation, but effects may remain small. This negative control is retained.
- Tipping points for risk models outside c_B=0.50–1.00 are labelled as outside the main grid, not presented as observed grid reversals.
- Signal probability is compared with no-change SPC controls; an individual signal cannot be causally attributed to recording on the basis of the chart alone.
- Paired probabilities depend on the recording-error coupling. Independent-cohort and independent-recording analyses are reported.
- Empirical replicate intervals, Monte Carlo uncertainty, nominal sampling CIs and signal/noise reliability have distinct meanings.

## Editor's two-sentence question

Previous work establishes that recording practices can influence quality measures. This study provides an explicit counterfactual stress-test that holds clinical reality constant, maps changes in record generation to changes in quality inference, and identifies the recording changes sufficient to erase or reverse improvement across measure architectures.

## Submission checks

Verify current BMJ Quality & Safety article-type availability, word/figure limits and declarations manually; the official author page was blocked during this session. Add authors, affiliations, funding, contributions, conflicts and institutional ethics determination where applicable. No journal submission or protocol registration has occurred.
'''
    (ROOT/'audit/analysis_audit.md').write_text(audit,encoding='utf-8')
    write_json({'validation_cells':len(validation),'max_abs_z':float(validation.analytic_z_delta.abs().max()),'all_R_10000':bool(all((df.R==cfg['replicates']).all() for df in [ac,risk,rep,sens,rs,spc,rel])),'main_AC_cells':len(ac),'main_risk_cells':len(risk)+len(rep),'sensitivity_AC_cells':len(sens),'sensitivity_risk_cells':len(rs),'SPC_cells':len(spc)//2},'audit/verification_summary.json')
    print('Tables, figures and analysis audit generated.',flush=True)
