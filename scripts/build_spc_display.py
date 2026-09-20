"""Display the complete, already executed SPC grid; no new simulations."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from src.common import ROOT,config
from src.reporting import style,save,md_table,BLUE,RED

def build_spc_display():
    cfg=config();n=cfg['spc']['monthly_N']
    data=pd.read_csv(ROOT/'results/main/spc.csv')
    main=data[data.n==n].copy();rows=[]
    endpoints=cfg['spc']['endpoints']
    for kind in ['abrupt','gradual']:
        for end in endpoints:
            selected=main[(main.kind==kind)&np.isclose(main.end,end)].set_index('rule')
            assert len(selected)==2
            primary=selected.loc['eight_same_side'];secondary=selected.loc['three_sigma']
            rows.append({'Change pattern':kind,'Final sensitivity (%)':end*100,
                         'Eight-point signal (%)':100*primary.signal_probability,
                         'Primary MCSE (pp)':100*primary.mcse,
                         'Three-sigma signal (%)':100*secondary.signal_probability,
                         'Secondary MCSE (pp)':100*secondary.mcse})
    table=pd.DataFrame(rows)
    table.to_csv(ROOT/'tables/S9_SPC_dose_response.csv',index=False,float_format='%.6g')
    md_table(table,'tables/S9_SPC_dose_response.md')
    caption=('All probabilities refer to at least one signal during the 24 follow-up months (months 25–48), '
             'with limits estimated from months 1–24 and frozen. Runs restart at month 25. '
             'These are monitoring-horizon probabilities, not single-look error probabilities. '
             f'Monthly n={n:,}; R=10,000 per scenario. The 80% endpoint is a no-change control. '
             'Abrupt and gradual controls use independent streams; their small difference is Monte Carlo variation. '
             'No new simulations were run to construct this table or figure.')
    path=ROOT/'tables/S9_SPC_dose_response.md'
    path.write_text('# Complete main SPC response across recording endpoints\n\n'+path.read_text(encoding='utf-8')+'\n'+caption+'\n',encoding='utf-8')
    style();fig,axes=plt.subplots(1,2,figsize=(10,4.6),sharey=True)
    for ax,rule,title in zip(axes,['eight_same_side','three_sigma'],['A  Eight consecutive same-side points','B  One point beyond three-sigma limits']):
        for kind,color,marker in [('abrupt',BLUE,'o'),('gradual',RED,'s')]:
            d=main[(main.rule==rule)&(main.kind==kind)].sort_values('end')
            p=d.signal_probability*100
            ax.errorbar((d.end-cfg['safety']['se_A'])*100,p,
                        yerr=[p-d.wilson_lo*100,d.wilson_hi*100-p],
                        color=color,marker=marker,capsize=3,lw=1.7,label=kind.capitalize())
        ax.set_title(title,loc='left',fontsize=10);ax.set_xlabel('Final ascertainment change (percentage points)')
        ax.set_xticks((np.array(endpoints)-cfg['safety']['se_A'])*100)
        ax.set_ylim(0,100);ax.grid(axis='y',alpha=.15);ax.legend(fontsize=9)
    axes[0].set_ylabel('At least one signal during 24-month follow-up (%)')
    fig.text(.5,.015,'Clinical risk remains 5%; baseline ascertainment is 80%. Zero change is the no-change comparator.\nError bars: Wilson 95% Monte Carlo intervals. All prespecified main endpoints are displayed.',ha='center',fontsize=9)
    fig.tight_layout(rect=[0,.1,1,1]);save(fig,'figureS5_spc_dose_response')
    return path.read_text(encoding='utf-8')

if __name__=='__main__':build_spc_display()
