"""GSE183320 壁细胞亚聚类、状态 panel 与供体级 DPT 趋势。

先运行 prepare_expression.py，再运行本文件。输出不采用细胞作为独立统计样本。
只有作者 log-normalized 矩阵，因此不运行要求 raw counts 的 DESeq2/edgeR。
"""
from pathlib import Path
import os
os.environ.setdefault('NUMBA_NUM_THREADS','4')
os.environ.setdefault('OMP_NUM_THREADS','4')
import json, importlib.metadata as im, warnings
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse, stats
from scipy.sparse.csgraph import connected_components
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

HERE=Path(__file__).resolve().parent
INPUT=HERE/'GSE183320_mural_EC_macrophage_normalized.h5ad'
OUT=HERE/'results'
SEED=20261013
PANELS={
 'contractile':['ACTA2','MYH11','TAGLN','CNN1','SMTN','MYL9','CALD1','LMOD1','ACTG2'],
 'matrix_fibrosis':['COL1A1','COL1A2','COL3A1','FN1','POSTN','CTGF','CCN2','SERPINE1','TIMP1','LOX','COL6A1'],
 'SASP_associated':['IL6','CXCL8','CCL2','CXCL1','CXCL2','MMP3','MMP9','IGFBP3','IGFBP7','GDF15'],
 'inflammatory':['IL1B','TNF','NFKBIA','RELA','ICAM1','VCAM1','CXCL10','CCL5','PTGS2','STAT1'],
}
# CTGF 与 CCN2 是同一基因的新旧 symbol：保留数据中实际存在的一种，不能重复加权。

def bh(p):
    p=np.asarray(p,float); q=np.full(p.shape,np.nan); ok=np.isfinite(p)
    v=p[ok]; order=np.argsort(v)
    if len(v):
        w=np.minimum.accumulate((v[order]*len(v)/np.arange(1,len(v)+1))[::-1])[::-1]
        z=np.empty(len(v));z[order]=np.minimum(w,1);q[ok]=z
    return q

def dense_gene(a,g): return a[:,g].X.toarray().ravel() if sparse.issparse(a.X) else np.asarray(a[:,g].X).ravel()

def partial_rank(x,y,z=None):
    x=np.asarray(x,float);y=np.asarray(y,float)
    if np.ptp(x)==0 or np.ptp(y)==0: return np.nan
    rx=stats.rankdata(x);ry=stats.rankdata(y)
    if z is not None:
        rz=stats.rankdata(np.asarray(z,float)); D=np.column_stack([np.ones(len(x)),rz])
        rx=rx-D@np.linalg.lstsq(D,rx,rcond=None)[0]
        ry=ry-D@np.linalg.lstsq(D,ry,rcond=None)[0]
    if np.std(rx)<1e-10 or np.std(ry)<1e-10: return np.nan
    return float(np.corrcoef(rx,ry)[0,1])

def signflip(values):
    v=np.asarray(values,float);v=v[np.isfinite(v)]
    if len(v)<5: return np.nan
    rng=np.random.default_rng(SEED)
    null=(rng.choice([-1,1],(9999,len(v)))*v).mean(1)
    return float((1+np.count_nonzero(np.abs(null)>=abs(v.mean())))/10000)

def donor_trends(a,ctype):
    rows=[];inclusion=[]
    for donor,g in a.obs.groupby('donor_id',observed=True):
        pos=a.obs_names.get_indexer(g.index)
        for predictor in ['dpt_primary','dpt_root_sensitivity']+['score_'+k for k in PANELS]:
            if predictor not in g: continue
            mask=np.isfinite(g[predictor].to_numpy(float)); gg=g.iloc[np.where(mask)[0]]
            eligible=len(gg)>=30 and (gg.PDE5A>0).sum()>=5
            inclusion.append({'celltype':ctype,'donor_id':str(donor),'predictor':predictor,
               'n_finite_cells':len(gg),'n_PDE5A_detected':int((gg.PDE5A>0).sum()),'eligible':eligible})
            if not eligible: continue
            for adjustment in ['unadjusted','nCount_adjusted']:
                rho=partial_rank(gg.PDE5A,gg[predictor],np.log1p(gg.nCount_RNA) if adjustment=='nCount_adjusted' else None)
                rows.append({'celltype':ctype,'donor_id':str(donor),'disease_status':str(gg.disease_status.iloc[0]),
                 'predictor':predictor,'adjustment':adjustment,'n_cells':len(gg),'rho':rho})
    pd.DataFrame(inclusion).to_csv(OUT/f'{ctype}_donor_trend_inclusion.csv',index=False)
    return rows

def marker_contrasts(a,ctype):
    # Gene family prescreen = >=10% detected in either group. 这里以整类 >=10%检出预筛，
    # 不看任何 p 值，完整所有检验都导出，之后统一 BH。
    genes=np.asarray(a.var_names)[np.asarray((a.X>0).mean(0)).ravel()>=.1]
    expr=a[:,genes].X
    rows=[]
    labels=a.obs['leiden_1.5'].astype(str)
    for cl in sorted(labels.unique(),key=int):
        deltas=[]; detdeltas=[]; donors=[]
        for donor in a.obs.donor_id.unique():
            m=(a.obs.donor_id==donor).to_numpy(); hi=m&(labels==cl).to_numpy();lo=m&(labels!=cl).to_numpy()
            if hi.sum()<10 or lo.sum()<20:continue
            deltas.append(np.asarray(expr[hi].mean(0)-expr[lo].mean(0)).ravel())
            detdeltas.append(np.asarray((expr[hi]>0).mean(0)-(expr[lo]>0).mean(0)).ravel())
            donors.append(str(donor))
        if len(donors)<5:continue
        D=np.asarray(deltas); E=np.asarray(detdeltas)
        # 零差全列不给伪 p=0；向量化 asymptotic Wilcoxon 并保留正负方向。
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            p=stats.wilcoxon(D,axis=0,zero_method='wilcox',alternative='two-sided',method='asymptotic').pvalue
        p[np.all(D==0,axis=0)]=1
        for j,gene in enumerate(genes):
            rows.append({'celltype':ctype,'cluster':cl,'gene':gene,'n_donors':len(donors),
             'mean_log_normalized_difference':D[:,j].mean(),'mean_detection_difference':E[:,j].mean(),
             'fraction_donors_positive_difference':(D[:,j]>0).mean(),'p_wilcoxon':p[j]})
    return pd.DataFrame(rows,columns=['celltype','cluster','gene','n_donors','mean_log_normalized_difference',
           'mean_detection_difference','fraction_donors_positive_difference','p_wilcoxon'])

def run_type(full,ctype):
    a=full[full.obs.celltype==ctype].copy()
    a.obs['donor_id']=a.obs.donor_id.astype(str)
    a.obs['PDE5A']=dense_gene(a,'PDE5A')
    # 保留作者元数据，用于每 donor 数量、技术深度审计。
    coverage=[]
    for name,gl in PANELS.items():
        present=[g for g in gl if g in a.var_names]
        if 'CTGF' in present and 'CCN2' in present: present.remove('CTGF')
        coverage.append({'celltype':ctype,'panel':name,'requested':'|'.join(gl),'present':'|'.join(present),
                         'n_present':len(present),'n_detected':sum((a[:,g].X>0).sum()>0 for g in present)})
        if len(present)<5: raise RuntimeError(f'Panel {name} has fewer than 5 measured genes')
        # 不让 PDE5A 或其它 panel 成员充当本模块的control，避免评分自关联。
        # 当前panel成员保留在gene_pool中供分箱，Scanpy会从control集合删除gene_list。
        forbidden=set(['PDE5A']+sum(PANELS.values(),[]))-set(present)
        pool=[gene for gene in a.var_names if gene not in forbidden]
        sc.tl.score_genes(a,present,score_name='score_'+name,ctrl_size=50,n_bins=25,
                         gene_pool=pool,random_state=SEED,use_raw=False,ctrl_as_ref=False)
    pd.DataFrame(coverage).to_csv(OUT/f'{ctype}_panel_coverage.csv',index=False)
    # graph genes 排除所有待检验基因/模块及技术基因；表达值只在 graph copy 中scale。
    sc.pp.filter_genes(a,min_cells=10)
    exclude=set(['PDE5A']+sum(PANELS.values(),[]))
    g=a[:,[v not in exclude and not v.startswith(('MT-','RPL','RPS')) for v in a.var_names]].copy()
    sc.pp.highly_variable_genes(g,flavor='seurat',n_top_genes=min(2000,g.n_vars),batch_key='donor_id',subset=False)
    g.var.loc[g.var.highly_variable].to_csv(OUT/f'{ctype}_HVG.csv')
    graph=g[:,g.var.highly_variable].copy()
    sc.pp.scale(graph,max_value=10)
    sc.tl.pca(graph,n_comps=min(30,graph.n_obs-1,graph.n_vars-1),svd_solver='arpack',random_state=SEED)
    sc.pp.neighbors(graph,n_neighbors=15,n_pcs=30,random_state=SEED)
    sc.tl.umap(graph,random_state=SEED)
    for res in [.8,1.2,1.5,2.0]:
        key=f'leiden_{res}'
        sc.tl.leiden(graph,resolution=res,random_state=SEED,key_added=key,flavor='leidenalg',directed=False)
        a.obs[key]=graph.obs[key]
    a.obsm['X_pca']=graph.obsm['X_pca'];a.obsm['X_umap']=graph.obsm['X_umap']
    a.obsp['connectivities']=graph.obsp['connectivities'];a.obsp['distances']=graph.obsp['distances']
    a.uns['neighbors']=graph.uns['neighbors']
    a.obs['UMAP_sub1']=a.obsm['X_umap'][:,0];a.obs['UMAP_sub2']=a.obsm['X_umap'][:,1]
    composition=pd.crosstab(a.obs.donor_id,a.obs['leiden_1.5'])
    composition.to_csv(OUT/f'{ctype}_cluster_donor_counts.csv')
    support=[]
    for cl in a.obs['leiden_1.5'].cat.categories:
        obs=a.obs.loc[a.obs['leiden_1.5']==cl]; counts=composition[cl]
        support.append({'celltype':ctype,'cluster':cl,'n_cells':len(obs),'n_donors_any':(counts>0).sum(),
         'n_donors_ge10':(counts>=10).sum(),'max_donor_fraction':counts.max()/counts.sum(),
         'donor_supported':len(obs)>=50 and (counts>=10).sum()>=5 and counts.max()/counts.sum()<=.4,
         'pooled_PDE5A_positive_pct':100*(obs.PDE5A>0).mean()})
    pd.DataFrame(support).to_csv(OUT/f'{ctype}_cluster_support.csv',index=False)
    # cell×donor 的依赖不能通过高分辨率解决：直接输出 donor mixing 描述。
    stability=[]
    for res in [.8,1.2,1.5,2.0]:
        label=a.obs[f'leiden_{res}'].astype(str)
        stability.append({'celltype':ctype,'resolution':res,'n_clusters':label.nunique(),
          'ARI_vs_primary':adjusted_rand_score(a.obs['leiden_1.5'],label),
          'NMI_with_donor':normalized_mutual_info_score(a.obs.donor_id,label)})
    pd.DataFrame(stability).to_csv(OUT/f'{ctype}_resolution_sensitivity.csv',index=False)
    # DPT 图只看转录邻接，不根据PDE5A选根/选方向。
    ncomp,comp=connected_components(graph.obsp['connectivities'],directed=False)
    big=np.bincount(comp).argmax(); keep=comp==big
    audit={'celltype':ctype,'n_cells':a.n_obs,'n_donors':a.obs.donor_id.nunique(),
           'n_graph_components':ncomp,'largest_component_fraction':float(keep.mean()),'seed':SEED}
    a.obs['dpt_primary']=np.nan;a.obs['dpt_root_sensitivity']=np.nan
    if keep.mean()>=.8:
        q=graph[keep].copy()
        q.obs['score_contractile']=a.obs.loc[q.obs_names,'score_contractile']
        # 最高10%中的PC质心邻近点，不使用PDE5A和病种标签。
        top=np.flatnonzero(q.obs.score_contractile>=q.obs.score_contractile.quantile(.9))
        center=q.obsm['X_pca'][top].mean(0)
        ordered=top[np.argsort(np.square(q.obsm['X_pca'][top]-center).sum(1))]
        sc.tl.diffmap(q,n_comps=15,random_state=SEED)
        roots=[]
        for i,key in enumerate(['dpt_primary','dpt_root_sensitivity']):
            q.uns['iroot']=int(ordered[i]);sc.tl.dpt(q,n_dcs=10)
            a.obs.loc[q.obs_names,key]=q.obs.dpt_pseudotime.to_numpy()
            roots.append(str(q.obs_names[ordered[i]]))
        audit['root_cells']=roots
        audit['root_sensitivity_spearman']=float(stats.spearmanr(a.obs.dpt_primary,a.obs.dpt_root_sensitivity,nan_policy='omit').statistic)
        audit['trajectory_label']='Contractile-reference transcriptional distance; not aging time'
    else:audit['dpt_status']='STOP: largest component below 80%'
    a.obs.to_csv(OUT/f'{ctype}_cell_annotations.csv.gz')
    a.write_h5ad(OUT/f'{ctype}_subclustered.h5ad',compression='gzip')
    (OUT/f'{ctype}_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    # 亚群里的PDE5A/模块只按 donor汇总，避免cell-pooled均值误读。
    d=a.obs.copy();d['positive']=d.PDE5A>0
    summaries=d.groupby(['donor_id','disease_status','leiden_1.5'],observed=True).agg(
     n_cells=('PDE5A','size'),PDE5A_mean=('PDE5A','mean'),PDE5A_pct=('positive',lambda x:100*x.mean()),
     **{k:(k,'mean') for k in ['score_'+name for name in PANELS]}).reset_index()
    summaries['eligible_ge10']=summaries.n_cells>=10
    summaries.to_csv(OUT/f'{ctype}_donor_cluster_summary.csv',index=False)
    marker_contrasts(a,ctype).to_csv(OUT/f'{ctype}_marker_contrasts_uncorrected.csv.gz',index=False)
    print(json.dumps(audit,indent=2),flush=True)
    return donor_trends(a,ctype)

def summarize_trends(df):
    rows=[];loo=[]
    for keys,g in df.groupby(['celltype','predictor','adjustment']):
        for cohort in ['all','control']:
            h=g if cohort=='all' else g[g.disease_status.astype(str).str.lower()=='control']
            v=h.rho.dropna().to_numpy(float);n=len(v)
            if not n:continue
            mean=v.mean();sd=v.std(ddof=1) if n>1 else np.nan
            width=stats.t.ppf(.975,n-1)*sd/np.sqrt(n) if n>1 else np.nan
            rows.append(dict(zip(['celltype','predictor','adjustment'],keys))|
             {'cohort':cohort,'n_donors':n,'mean_partial_rho':mean,'sd_donor_rho':sd,
             'mean_rho_ci95_lower':max(-1,mean-width),'mean_rho_ci95_upper':min(1,mean+width),
             'positive_donors':int((v>0).sum()),'p_signflip':signflip(v)})
            for _,x in h.iterrows():
                other=h.loc[h.donor_id!=x.donor_id,'rho'].dropna()
                if len(other):loo.append(dict(zip(['celltype','predictor','adjustment'],keys))|
                  {'cohort':cohort,'omitted_donor':x.donor_id,'n_remaining':len(other),'mean_rho':other.mean()})
    out=pd.DataFrame(rows)
    # 所有 2types×5 主要预测变量一起BH，不选显著模块。根敏感性单列，不并入主 family。
    out['q_BH']=np.nan
    for cohort in ['all','control']:
        for adj in ['unadjusted','nCount_adjusted']:
            m=(out.cohort==cohort)&(out.adjustment==adj)&(out.predictor!='dpt_root_sensitivity')
            out.loc[m,'q_BH']=bh(out.loc[m,'p_signflip'])
    out.to_csv(OUT/'donor_trend_summary.csv',index=False)
    pd.DataFrame(loo).to_csv(OUT/'donor_trend_leave_one_out.csv',index=False)

def finalize():
    frames=[]
    for ct in ['SMC','pericyte']:
        try: frames.append(pd.read_csv(OUT/f'{ct}_marker_contrasts_uncorrected.csv.gz'))
        except pd.errors.EmptyDataError:
            # 完全没有亚群通过固定最少供体门槛也是结果；不能为产出marker降低门槛。
            empty=pd.DataFrame(columns=['celltype','cluster','gene','n_donors','mean_log_normalized_difference',
               'mean_detection_difference','fraction_donors_positive_difference','p_wilcoxon'])
            empty.to_csv(OUT/f'{ct}_marker_contrasts_uncorrected.csv.gz',index=False);frames.append(empty)
    markers=pd.concat(frames,ignore_index=True)
    markers['q_BH_all_cluster_gene_tests']=bh(markers.p_wilcoxon)
    markers['descriptive_marker']=((markers.q_BH_all_cluster_gene_tests<.05)&
       (markers.mean_log_normalized_difference>0)&(markers.mean_detection_difference>=.1))
    markers.to_csv(OUT/'all_donor_paired_markers.csv.gz',index=False)
    markers[markers.descriptive_marker].sort_values(['celltype','cluster','q_BH_all_cluster_gene_tests']).to_csv(OUT/'candidate_marker_table.csv',index=False)
    versions={k:im.version(k) for k in ['scanpy','anndata','numpy','pandas','scipy','scikit-learn','igraph','leidenalg','umap-learn']}
    (OUT/'software_versions.json').write_text(json.dumps(versions,indent=2),encoding='utf-8')
    print('Completed donor-aware analysis',flush=True)

def main():
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--finalize-only',action='store_true')
    args=parser.parse_args();OUT.mkdir(exist_ok=True)
    if args.finalize_only:finalize();return
    full=ad.read_h5ad(INPUT)
    # 先切小矩阵释放EC/macrophage，SMC+pericyte共2666细胞。
    mural=full[full.obs.celltype.isin(['SMC','pericyte'])].copy();del full
    rows=[]
    for ctype in ['SMC','pericyte']:rows.extend(run_type(mural,ctype))
    trend=pd.DataFrame(rows);trend.to_csv(OUT/'donor_trends.csv',index=False);summarize_trends(trend)
    finalize()

if __name__=='__main__':main()
