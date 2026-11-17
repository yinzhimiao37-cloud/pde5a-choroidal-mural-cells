"""真实供体结果绘图；PDF/SVG可编辑文本，误差线为供体间样本SD。"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.transforms import ScaledTranslation

HERE=Path(__file__).resolve().parent
OUT=HERE/'figures';OUT.mkdir(exist_ok=True)
# Optional display-only alignment helper; not part of statistical computation.
try:
    from audit_panel_alignment import require_matplotlib_panel_alignment
except ModuleNotFoundError:
    def require_matplotlib_panel_alignment(*args, **kwargs):
        import warnings
        warnings.warn("Optional panel-alignment helper unavailable; alignment audit skipped. Inspect the exported figure manually.", RuntimeWarning)


plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
    'font.size':7.5,'axes.labelsize':8,'axes.titlesize':8.5,'xtick.labelsize':7.2,'ytick.labelsize':7.2,
    'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':0.7,
    'pdf.fonttype':42,'svg.fonttype':'none','legend.frameon':False,'savefig.facecolor':'white'})
BLUE='#36799E';TEAL='#529C90';ORANGE='#BD7642';GRAY='#9C9C9C'
COL={'SMC':BLUE,'pericyte':TEAL,'fibroblast':GRAY,'endothelial_combined':'#B6B6B6'}
LAB={'SMC':'SMC','pericyte':'Pericyte','fibroblast':'Fibroblast','endothelial_combined':'Endothelium'}
# 随机数仅用于水平防重叠jitter；所有纵坐标、样本量和误差线来自真实数据。
rng=np.random.default_rng(20261013)
by=pd.read_csv(HERE/'GSE183320/donor_celltype_PDE5A.csv')
th=pd.read_csv(HERE/'GSE183320/threshold_cohort_summary.csv')
loo=pd.read_csv(HERE/'GSE183320/leave_one_donor_out.csv')
md=pd.read_csv(HERE/'GSE183320/depth_matched_donor_summary.csv')

def panel(ax,label,title):
    ax.set_title(title,loc='left',pad=9,fontweight='normal')
    ax.text(0,1,label,transform=ax.transAxes+ScaledTranslation(-21/72,10/72,ax.figure.dpi_scale_trans),
            fontsize=9,fontweight='bold',va='bottom',ha='left')
    ax.tick_params(length=3,width=0.7)

def save(fig,name):
    fig.canvas.draw()
    require_matplotlib_panel_alignment(fig,json_out=OUT/(name+'.alignment.json'),
        overlay_svg=OUT/(name+'.alignment.svg'),tolerance_pt=1.5,gutter_tolerance_pt=1.5,
        require_panel_labels=True,strict=True)
    stem=OUT/name
    fig.savefig(f'{stem}.pdf')
    fig.savefig(f'{stem}.svg')
    fig.savefig(f'{stem}.png',dpi=300)
    fig.savefig(f'{stem}.tiff',dpi=600)
    plt.close(fig)

def dotgroups(ax,frame,types,ycol='pct_positive'):
    labels=[]
    for i,t in enumerate(types):
        q=frame[frame.celltype.eq(t)].sort_values('donor_id')
        vals=q[ycol].to_numpy()
        ax.scatter(i+rng.uniform(-0.15,0.08,len(vals)),vals,s=20,c=COL[t],alpha=.85,
                   edgecolors='white',linewidths=.35,zorder=3)
        if len(vals):ax.errorbar(i+.21,vals.mean(),yerr=vals.std(ddof=1),fmt='o',c='#202020',ms=3.1,
                               capsize=3,elinewidth=1,zorder=4)
        labels.append(f'{LAB[t]}\nn = {len(vals)}')
    ax.set_xticks(range(len(types)),labels)
    ax.set_xlim(-.45,len(types)-.45)
    ax.set_ylim(-3,103);ax.set_yticks([0,25,50,75,100])
    ax.set_ylabel('PDE5A-positive cells (%)')
    ax.grid(axis='y',color='#E9E9E9',lw=.5,zorder=0)

# 定位、健康限制、来源复用和深度平衡各回答一个不同问题。
fig,axes=plt.subplots(2,2,figsize=(7.086614,6.220472))
fig.subplots_adjust(left=.10,right=.98,bottom=.10,top=.88,wspace=.34,hspace=.52)
primary=by[by.cells_n.ge(20)]
types=['SMC','pericyte','fibroblast','endothelial_combined']
for ax,label,title,sub in [
    (axes[0,0],'a','All donors',primary),
    (axes[0,1],'b','Control donors',primary[primary.disease_status.eq('control')]),
    (axes[1,0],'c','Known GSE135922 overlap excluded',primary[~primary.donor_id.isin([1,2,11,12])])]:
    panel(ax,label,title);dotgroups(ax,sub,types)
ax=axes[1,1];panel(ax,'d','After donor-wise depth-bin matching')
for i,t in enumerate(['SMC','pericyte']):
    q=md[md.celltype.eq(t)];v=q.difference_pp
    ax.scatter(i+rng.uniform(-.13,.08,len(q)),v,s=22,c=COL[t],edgecolors='white',linewidths=.3,zorder=3)
    ax.errorbar(i+.2,v.mean(),yerr=v.std(ddof=1),fmt='o',c='#222222',ms=3.2,capsize=3,elinewidth=1)
ax.set_xticks([0,1],[f'SMC − endothelium\nn = {sum(md.celltype.eq("SMC"))}',
                    f'Pericyte − endothelium\nn = {sum(md.celltype.eq("pericyte"))}'])
ax.set_ylabel('Difference in positivity (percentage points)')
ax.set_xlim(-.4,1.5);ax.set_ylim(-3,103);ax.set_yticks([0,25,50,75,100]);ax.axhline(0,color='#777777',lw=.7)
ax.grid(axis='y',color='#E9E9E9',lw=.5)
fig.text(.10,.975,'GSE183320  |  PDE5A cellular localization',fontsize=10,weight='bold',va='top')
save(fig,'Figure1_GSE183320_donor_localization')

# 阈值、单供体影响与技术相关性，保持完整的敏感性结果。
fig,axes=plt.subplots(2,2,figsize=(7.086614,6.417323))
fig.subplots_adjust(left=.10,right=.98,bottom=.11,top=.90,wspace=.36,hspace=.70)
levels=[1,5,10,20,30,50]
for ax,label,typ in [(axes[0,0],'a','SMC'),(axes[0,1],'b','pericyte')]:
    panel(ax,label,LAB[typ]+' cell-count sensitivity')
    ns=[]
    for cohort,color,offset,display in [('all',BLUE,-.06,'All donors'),('control',ORANGE,.06,'Control')]:
        q=th[th.celltype.eq(typ)&th.cohort.eq(cohort)].set_index('min_cells').loc[levels]
        ax.errorbar(np.arange(6)+offset,q.mean_percent,yerr=q.sample_sd_percent,fmt='o-',
                    color=color,ms=3.0,lw=1,capsize=2,label=display)
        ns.append(q.donors_n.to_list())
    ax.set_xticks(range(6),[f'{n}\n{a}/{b}' for n,a,b in zip(levels,*ns)])
    ax.set_xlabel('Minimum cells per donor\nDonors: all / control')
    ax.set_ylabel('PDE5A-positive cells (%)');ax.set_ylim(0,110);ax.set_yticks([0,25,50,75,100])
    ax.grid(axis='y',color='#E9E9E9',lw=.5)
handles=[Line2D([0],[0],color=BLUE,marker='o',ms=3,label='All donors'),
         Line2D([0],[0],color=ORANGE,marker='o',ms=3,label='Control')]
fig.legend(handles=handles,loc='upper right',bbox_to_anchor=(.98,.973),ncol=2,fontsize=7.5)
ax=axes[1,0];panel(ax,'c','Leave-one-donor-out influence')
for i,(typ,cohort) in enumerate([('SMC','all'),('SMC','control'),('pericyte','all'),('pericyte','control')]):
    q=loo[loo.celltype.eq(typ)&loo.cohort.eq(cohort)]
    ax.scatter(i+rng.uniform(-.12,.12,len(q)),q.shift_percentage_points,s=16,
               c=BLUE if cohort=='all' else ORANGE,edgecolors='white',linewidths=.3)
ax.axhline(0,c='#777777',lw=.7)
ax.set_xticks(range(4),['SMC\nAll, n = 17','SMC\nControl, n = 8','Pericyte\nAll, n = 20','Pericyte\nControl, n = 10'])
ax.set_ylabel('Change in mean positivity (pp)');ax.set_ylim(-8,8)
ax=axes[1,1];panel(ax,'d','Detection-depth association in SMC')
q=primary[primary.celltype.eq('SMC')]
ax.scatter(q.median_nFeature_RNA,q.pct_positive,s=25,c=BLUE,edgecolors='white',linewidths=.4)
ax.set_xlabel('Median detected genes per donor')
ax.set_ylabel('PDE5A-positive SMC (%)');ax.set_ylim(0,110)
ax.text(.03,.96,'Spearman ρ = 0.578\nPermutation P = 0.0168; BH q = 0.0336\nn = 17 donors',
        transform=ax.transAxes,va='top',fontsize=7)
fig.text(.10,.98,'GSE183320  |  Robustness and technical sensitivity',fontsize=10,weight='bold',va='top')
save(fig,'Figure2_GSE183320_sensitivity')

# GSE135922区域是供体内配对；不能把14个区域称14位供体。
g=HERE/'GSE135922'
d=pd.read_csv(g/'PDE5A_mural_donor_pooled_regions.csv')
r=pd.read_csv(g/'PDE5A_mural_donor_region.csv')
fig,axes=plt.subplots(1,2,figsize=(7.086614,3.582677))
fig.subplots_adjust(left=.10,right=.98,bottom=.20,top=.76,wspace=.35)
ax=axes[0];panel(ax,'a','Donor-level mural-like localization')
groups=[d,d[d.donor.isin([1,2,3])]]
for i,q in enumerate(groups):
    for j,row in enumerate(q.itertuples()):
        ax.scatter(i+(j-(len(q)-1)/2)*.055,row.positivity_percent,s=30,
                   c=BLUE if row.experiment=='unenriched' else ORANGE,edgecolors='white',linewidths=.4)
    v=q.positivity_percent
    ax.errorbar(i+.27,v.mean(),yerr=v.std(ddof=1),fmt='o',ms=3,c='#222222',capsize=3,elinewidth=1)
ax.set_xticks([0,1],['All donors\nn = 7','Known reused donors excluded\nn = 3'])
ax.set_ylabel('PDE5A-positive mural-like cells (%)');ax.set_ylim(0,65);ax.set_xlim(-.4,1.55)
ax=axes[1];panel(ax,'b','Within-donor regional comparison')
for donor,q in r.groupby('donor'):
    q=q.set_index('region').loc[['macula','peripheral']]
    eligible=q.cells_n.ge(20).all()
    col=BLUE if q.experiment.iloc[0]=='unenriched' else ORANGE
    ax.plot([0,1],q.positivity_percent,color=col,lw=.9,alpha=.75,ls='-' if eligible else '--')
    for x,(_,row) in enumerate(q.iterrows()):
        ax.scatter(x,row.positivity_percent,s=29,facecolors=col if row.cells_n>=20 else 'white',
                   edgecolors=col,linewidths=.9,zorder=3)
ax.set_xticks([0,1],['Macula','Peripheral']);ax.set_xlim(-.25,1.25);ax.set_ylim(0,65)
ax.set_ylabel('PDE5A-positive mural-like cells (%)')
fig.text(.10,.96,'GSE135922  |  Donor and regional consistency',fontsize=10,weight='bold',va='top')
fig.legend(handles=[Line2D([0],[0],marker='o',color=BLUE,label='Unenriched, donors 1–3',ms=3),
                    Line2D([0],[0],marker='o',color=ORANGE,label='CD31-enriched, donors 4–7',ms=3)],
           loc='upper left',bbox_to_anchor=(.08,.895),ncol=2,fontsize=7.3)
fig.text(.10,.035,'Solid: 4 eligible pairs; dashed: ≥1 region <20 cells; hollow: region <20 cells.',fontsize=7,ha='left')
save(fig,'Figure3_GSE135922_donor_region')
print('Exported 3 figures with alignment checks.')
