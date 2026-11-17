"""绘制真实Xenium供体点和距离曲线；置换包络不是供体置信区间。"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.transforms import ScaledTranslation

H=Path(__file__).resolve().parent
O=H/'figures'
# Optional display-only alignment helper; not part of statistical computation.
try:
    from audit_panel_alignment import require_matplotlib_panel_alignment
except ModuleNotFoundError:
    def require_matplotlib_panel_alignment(*args, **kwargs):
        import warnings
        warnings.warn("Optional panel-alignment helper unavailable; alignment audit skipped. Inspect the exported figure manually.", RuntimeWarning)

plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
    'font.size':7.5,'axes.labelsize':8,'axes.titlesize':8.5,'xtick.labelsize':7.2,'ytick.labelsize':7.2,
    'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,
    'pdf.fonttype':42,'svg.fonttype':'none','legend.frameon':False,'savefig.facecolor':'white'})
don=pd.read_csv(H/'spatial_GSE318037/donor_spatial_summary.csv')
co=pd.read_csv(H/'spatial_GSE318037/cohort_spatial_summary.csv')
cv=pd.read_csv(H/'spatial_GSE318037/distance_gradient_source.csv')
def primary(x):return x[x.PDE5A_min_transcripts.eq(1)&x.scheme.eq('depth_density_stratified')]
don=primary(don);co=primary(co[co.cohort.eq('all')]);cv=primary(cv[cv.cohort.eq('all')])
fig,axs=plt.subplots(2,2,figsize=(7.086614,6.30))
fig.subplots_adjust(left=.105,right=.97,top=.87,bottom=.18,wspace=.35,hspace=.58)
for col,(typ,color) in enumerate([('SMC','#36799E'),('pericyte','#529C90')]):
    label='SMC' if typ=='SMC' else 'Pericyte'
    ax=axs[0,col];q=don[don.celltype.eq(typ)].sort_values('donor')
    r=co[co.celltype.eq(typ)].iloc[0]
    # 确定性水平偏移仅防遮挡；纵坐标全部来自真实供体结果。
    offsets=np.linspace(-.11,.11,len(q))
    for dx,(_,d) in zip(offsets,q.iterrows()):
        ys=[d.EC_neighbors_25um,d.EC_neighbors_25um_null_mean]
        ax.plot([dx,1+dx],ys,color='#C4C4C4',lw=.65,zorder=1)
        marker='o' if d.condition=='Control' else '^'
        ax.scatter(dx,ys[0],s=25,marker=marker,color=color,edgecolor='white',lw=.35,zorder=3)
        ax.scatter(1+dx,ys[1],s=25,marker=marker,color='#929292',edgecolor='white',lw=.35,zorder=3)
    for x,name in enumerate(['EC_neighbors_25um','EC_neighbors_25um_null_mean']):
        v=q[name]
        ax.errorbar(x+.22,v.mean(),yerr=v.std(ddof=1),fmt='o',c='#222222',ms=3,capsize=3,elinewidth=1,zorder=4)
    ax.set_xticks([0,1],['Observed','Conditional null'])
    ax.set_xlim(-.35,1.4);ax.set_ylim(0,5.2);ax.set_yticks([0,1,2,3,4])
    ax.set_ylabel('Endothelial neighbours within 25 µm')
    ax.text(.03,.97,f'n = {len(q)} donors; {int(r.positive_cells_n)} positive cells\nNES = {r.NES:.2f}; P = {r.permutation_p_two_sided:.4f}\nBH q = {r.BH_q_primary_family:.4f}',transform=ax.transAxes,va='top',fontsize=7.2)
    ax.set_title(label+' | fixed spatial endpoint',loc='left',pad=9)
    ax.grid(axis='y',color='#ECECEC',lw=.5,zorder=0)
    ax=axs[1,col];c=cv[cv.celltype.eq(typ)].sort_values('radius_um')
    # 彩色带：供体间SD；灰带：固定组织条件下的95%置换包络。
    x=c.radius_um.to_numpy();v=c.observed_fraction.to_numpy();sd=c.sd_across_donors.to_numpy()
    ax.fill_between(x,100*np.clip(v-sd,0,1),100*np.clip(v+sd,0,1),color=color,alpha=.13,lw=0)
    ax.fill_between(x,100*c.null_lower_025,100*c.null_upper_975,color='#777777',alpha=.18,lw=0)
    ax.plot(x,100*c.null_mean_fraction,color='#777777',ls='--',lw=1.1)
    ax.plot(x,100*v,color=color,lw=1.4,marker='o',ms=2.5)
    ax.set_xlim(0,100);ax.set_ylim(0,103);ax.set_xticks([0,25,50,75,100]);ax.set_yticks([0,25,50,75,100])
    ax.set_xlabel('Distance to nearest endothelial centroid (µm)')
    ax.set_ylabel('Positive mural cells within distance (%)')
    ax.set_title(label+' | distance distribution',loc='left',pad=9)
    ax.grid(axis='y',color='#ECECEC',lw=.5,zorder=0)
for ax,letter in zip(axs.flat,'abcd'):
    ax.text(0,1,letter,transform=ax.transAxes+ScaledTranslation(-21/72,10/72,fig.dpi_scale_trans),
        fontsize=9,fontweight='bold',va='bottom',ha='left')
    ax.tick_params(length=3,width=.7)
fig.text(.105,.98,'GSE318037 | PDE5A mural-cell proximity to endothelium',fontsize=10,weight='bold',va='top')
fig.legend(handles=[Line2D([],[],marker='o',ls='',color='#666666',ms=4,label='Non-CNV control tissue'),
    Line2D([],[],marker='^',ls='',color='#666666',ms=4,label='CNV donor')],loc='upper left',bbox_to_anchor=(.095,.945),ncol=2,fontsize=7.2)
fig.text(.105,.04,'a,b: paired donor points; black = mean ± donor SD. c,d: colour = observed mean ± donor SD;\ngrey = conditional null mean and 95% permutation envelope. Distances refer to cell centroids.',fontsize=7.1,va='bottom',linespacing=1.5)
fig.canvas.draw()
name='Figure4_GSE318037_spatial_context'
require_matplotlib_panel_alignment(fig,json_out=O/(name+'.alignment.json'),overlay_svg=O/(name+'.alignment.svg'),
    tolerance_pt=1.5,gutter_tolerance_pt=1.5,require_panel_labels=True,strict=True)
fig.savefig(O/(name+'.pdf'))
fig.savefig(O/(name+'.svg'))
fig.savefig(O/(name+'.png'),dpi=300)
fig.savefig(O/(name+'.tiff'),dpi=600)
plt.close(fig)
print('Saved',name)
