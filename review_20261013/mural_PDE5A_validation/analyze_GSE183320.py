"""事后供体级PDE5A壁细胞分析。阈值固定于protocol_locked.md。"""
from pathlib import Path
import csv, gzip, hashlib, json, platform
from datetime import datetime, timezone
import numpy as np
import pandas as pd
import scipy
from scipy.stats import rankdata

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
OUT = HERE/'GSE183320'
OUT.mkdir(exist_ok=True)
SEED = 20261013
THRESHOLDS = [1,5,10,20,30,50]
TYPES = ['SMC','pericyte','mural_combined','fibroblast','endothelial_combined']
EC = ['artery','vein','choriocapillaris','ec_unknown']
OVERLAP = {1,2,11,12}
PRIOR = {1,2,3,5,11,12}

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()

protocol=HERE/'protocol_locked.md'
audit={'started_utc':datetime.now(timezone.utc).isoformat(), 'protocol_sha256':sha(protocol),
       'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,
       'scipy':scipy.__version__, 'seed':SEED,'post_hoc':True,
       'positivity':'published normalized PDE5A > 0; not raw counts'}
meta_path=PROJECT/'data/raw/GSE183320_all_metadata.csv.gz'
expr_path=PROJECT/'data/processed/sc_choroid_focus_matrix.csv'
raw_path=PROJECT/'data/raw/GSE183320_all_normalized_expression_data.csv.gz'
meta=pd.read_csv(meta_path)
expr=pd.read_csv(expr_path,usecols=['cell_barcode','PDE5A'])
meta['cell_barcode']=meta.barcode.str.replace('-','.',regex=False)
assert meta.cell_barcode.is_unique and expr.cell_barcode.is_unique
dat=meta.merge(expr,on='cell_barcode',validate='one_to_one',how='outer',indicator=True)
assert len(dat)==30416 and dat['_merge'].eq('both').all()
dat=dat.drop(columns='_merge')
assert dat[['PDE5A','nCount_RNA','nFeature_RNA']].notna().all().all()
assert np.isfinite(dat[['PDE5A','nCount_RNA','nFeature_RNA']].to_numpy()).all()
assert dat.PDE5A.ge(0).all() and dat.nCount_RNA.gt(0).all()
dat['donor_id']=dat.library.str.extract(r'donor_(\d+)_')[0].astype(int)
assert dat.groupby('donor_id').library.nunique().max()==1
assert dat.groupby('donor_id').disease_status.nunique().max()==1
dat['positive']=dat.PDE5A.gt(0)

# 直接从原始发表矩阵抽出PDE5A一行，验证旧focus表无条码或数值错位。
with gzip.open(raw_path,'rt',newline='') as f:
    reader=csv.reader(f)
    header=next(reader)
    for line in reader:
        if line[0].strip('"')=='PDE5A':
            source=pd.Series(np.array(line[1:],dtype=float),index=header[1:]);break
    else:raise ValueError('PDE5A missing from published matrix')
source.index=source.index.str.replace('-','.',regex=False)
assert source.index.is_unique
v=dat.cell_barcode.map(source)
assert v.notna().all()
audit['max_abs_expression_difference']=float(np.max(np.abs(v-dat.PDE5A)))
assert np.allclose(v,dat.PDE5A,rtol=1e-12,atol=1e-12)
audit['source_validation_passed']=True
audit['input_sha256']={str(p.relative_to(PROJECT)):sha(p) for p in [meta_path,expr_path,raw_path]}
donors=dat[['donor_id','library','disease_status','age','sex','eye']].drop_duplicates().sort_values('donor_id')
assert donors.donor_id.is_unique
donors['known_GSE135922_reuse']=donors.donor_id.isin(OVERLAP)
donors['any_prior_series_reuse']=donors.donor_id.isin(PRIOR)
donors.to_csv(OUT/'donor_metadata.csv',index=False)
audit['libraries_n']=len(donors)
audit['disease_status_counts']=donors.disease_status.value_counts().to_dict()

def mask_type(t):
    if t=='mural_combined':return dat.celltype.isin(['SMC','pericyte'])
    if t=='endothelial_combined':return dat.celltype.isin(EC)
    return dat.celltype.eq(t)

rows=[]
for typ in TYPES:
    sub=dat[mask_type(typ)]
    for donor in donors.itertuples():
        x=sub[sub.donor_id.eq(donor.donor_id)]
        n=len(x); k=int(x.positive.sum())
        rows.append(dict(donor_id=donor.donor_id,library=donor.library,disease_status=donor.disease_status,
            celltype=typ,cells_n=n,positive_n=k,pct_positive=100*k/n if n else np.nan,
            mean_normalized_expression=x.PDE5A.mean(),median_nCount_RNA=x.nCount_RNA.median(),
            median_nFeature_RNA=x.nFeature_RNA.median(),known_GSE135922_reuse=donor.known_GSE135922_reuse,
            any_prior_series_reuse=donor.any_prior_series_reuse,eligible_primary=n>=20))
by=pd.DataFrame(rows)
by.to_csv(OUT/'donor_celltype_PDE5A.csv',index=False)
dat[dat.celltype.isin(['SMC','pericyte','fibroblast']+EC)][['barcode','donor_id','library',
    'disease_status','celltype','PDE5A','positive','nCount_RNA','nFeature_RNA']].to_csv(OUT/'cell_source_data.csv.gz',index=False)

def cohort_mask(df,name):
    if name=='all':return pd.Series(True,index=df.index)
    if name=='control':return df.disease_status.eq('control')
    if name=='no_known_GSE135922_overlap':return ~df.donor_id.isin(OVERLAP)
    if name=='control_no_known_GSE135922_overlap':return df.disease_status.eq('control') & ~df.donor_id.isin(OVERLAP)
    if name=='new_series_only':return ~df.donor_id.isin(PRIOR)
    if name=='control_new_series_only':return df.disease_status.eq('control') & ~df.donor_id.isin(PRIOR)
    raise ValueError(name)

COHORTS=['all','control','no_known_GSE135922_overlap','control_no_known_GSE135922_overlap',
         'new_series_only','control_new_series_only']
summary=[];exclusion=[]
for cohort in COHORTS:
    base=by[cohort_mask(by,cohort)]
    for typ in TYPES:
        t=base[base.celltype.eq(typ)]
        for threshold in THRESHOLDS:
            q=t[t.cells_n.ge(threshold)]
            summary.append(dict(cohort=cohort,celltype=typ,min_cells=threshold,donors_n=len(q),
                mean_percent=q.pct_positive.mean(),sample_sd_percent=q.pct_positive.std(ddof=1),
                cells_n=int(q.cells_n.sum()),positive_n=int(q.positive_n.sum()),
                pooled_percent=100*q.positive_n.sum()/q.cells_n.sum() if q.cells_n.sum() else np.nan,
                median_cells_per_donor=q.cells_n.median(),min_cells_observed=q.cells_n.min(),
                max_cells_observed=q.cells_n.max(),eligible_donor_ids=';'.join(map(str,q.donor_id)),
                observed_donors_n=int(t.cells_n.gt(0).sum()),cohort_donors_n=len(t)))
            for r in t.itertuples():exclusion.append(dict(cohort=cohort,celltype=typ,min_cells=threshold,
                donor_id=r.donor_id,cells_n=r.cells_n,included=r.cells_n>=threshold,
                reason='included' if r.cells_n>=threshold else ('no_cells' if r.cells_n==0 else 'below_min_cells')))
summary=pd.DataFrame(summary)
summary.to_csv(OUT/'threshold_cohort_summary.csv',index=False)
pd.DataFrame(exclusion).to_csv(OUT/'inclusion_audit.csv',index=False)
primary=summary[summary.min_cells.eq(20)]
primary.to_csv(OUT/'primary_summary.csv',index=False)

# LOO只剔除一个完整供体，不把剩余细胞重新当独立重复。
loo=[]
for cohort in COHORTS:
    for typ in TYPES:
        q=by[cohort_mask(by,cohort)&by.celltype.eq(typ)&by.cells_n.ge(20)]
        if len(q)<3:continue
        baseline=q.pct_positive.mean()
        for r in q.itertuples():
            z=q[q.donor_id.ne(r.donor_id)]
            loo.append(dict(cohort=cohort,celltype=typ,removed_donor=r.donor_id,remaining_donors_n=len(z),
                baseline_mean_percent=baseline,loo_mean_percent=z.pct_positive.mean(),
                loo_sd_percent=z.pct_positive.std(ddof=1),shift_percentage_points=z.pct_positive.mean()-baseline))
loo=pd.DataFrame(loo);loo.to_csv(OUT/'leave_one_donor_out.csv',index=False)
ls=loo.groupby(['cohort','celltype']).agg(loo_mean_min=('loo_mean_percent','min'),
    loo_mean_max=('loo_mean_percent','max'),max_abs_shift_pp=('shift_percentage_points',lambda x:abs(x).max())).reset_index()
ls.to_csv(OUT/'leave_one_donor_out_summary.csv',index=False)

# 小n供体Spearman置换：在供体之间置换，不在细胞之间置换。
rng=np.random.default_rng(SEED)
cor=[]
for cohort in ['all','control']:
    for typ in ['SMC','pericyte']:
        q=by[cohort_mask(by,cohort)&by.celltype.eq(typ)&by.cells_n.ge(20)]
        for metric in ['median_nCount_RNA','median_nFeature_RNA']:
            n=len(q);rho=np.nan;p=np.nan
            x=rankdata(q[metric]);y=rankdata(q.pct_positive)
            if n>=4 and np.std(x)>0 and np.std(y)>0:
                x=(x-x.mean())/np.linalg.norm(x-x.mean())
                y=(y-y.mean())/np.linalg.norm(y-y.mean())
                rho=float(x@y)
                null=np.array([x@rng.permutation(y) for _ in range(9999)])
                p=float((1+np.count_nonzero(np.abs(null)>=abs(rho)-1e-12))/10000)
            cor.append(dict(cohort=cohort,celltype=typ,metric=metric,donors_n=n,spearman_rho=rho,
                permutation_p_two_sided=p,permutations=9999,seed=SEED))
cor=pd.DataFrame(cor)
cor['BH_q']=np.nan
for _,g in cor.groupby(['cohort','celltype']):
    good=g.permutation_p_two_sided.dropna().sort_values()
    if len(good):
        adj=np.minimum.accumulate((good.to_numpy()*len(good)/np.arange(1,len(good)+1))[::-1])[::-1]
        cor.loc[good.index,'BH_q']=np.minimum(adj,1)
cor.to_csv(OUT/'donor_depth_correlations.csv',index=False)

# 在供体内按深度五分位匹配细胞数量，不把归一化数值反算UMI。
matched=[];match_audit=[]
for donor in donors.itertuples():
    dd=dat[dat.donor_id.eq(donor.donor_id)]
    for typ in ['SMC','pericyte']:
        a=dd[dd.celltype.eq(typ)].copy();b=dd[dd.celltype.isin(EC)].copy()
        reason='eligible'
        if len(a)<20 or len(b)<20:reason='original_group_below_20'
        if reason!='eligible':
            match_audit.append(dict(donor_id=donor.donor_id,celltype=typ,mural_n=len(a),endothelial_n=len(b),matched_n=0,reason=reason));continue
        edges=np.unique(np.quantile(np.log1p(pd.concat([a.nCount_RNA,b.nCount_RNA])),np.linspace(0,1,6)))
        if len(edges)<3:reason='insufficient_depth_bins'
        else:
            edges[0]=-np.inf;edges[-1]=np.inf
            a['bin']=pd.cut(np.log1p(a.nCount_RNA),bins=edges,labels=False,include_lowest=True)
            b['bin']=pd.cut(np.log1p(b.nCount_RNA),bins=edges,labels=False,include_lowest=True)
            bins=[(a.index[a.bin.eq(j)].to_numpy(),b.index[b.bin.eq(j)].to_numpy()) for j in range(len(edges)-1)]
            nmatch=sum(min(len(aa),len(bb)) for aa,bb in bins)
            if nmatch<20:reason='matched_group_below_20'
        if reason!='eligible':
            match_audit.append(dict(donor_id=donor.donor_id,celltype=typ,mural_n=len(a),endothelial_n=len(b),matched_n=0,reason=reason));continue
        match_audit.append(dict(donor_id=donor.donor_id,celltype=typ,mural_n=len(a),endothelial_n=len(b),matched_n=nmatch,reason=reason))
        for rep in range(100):
            ia=[];ib=[]
            for aa,bb in bins:
                count=min(len(aa),len(bb))
                if count:
                    ia.extend(rng.choice(aa,count,replace=False));ib.extend(rng.choice(bb,count,replace=False))
            ap=100*dat.loc[ia,'positive'].mean();bp=100*dat.loc[ib,'positive'].mean()
            la=np.log1p(dat.loc[ia,'nCount_RNA'].to_numpy())
            lb=np.log1p(dat.loc[ib,'nCount_RNA'].to_numpy())
            pooled_sd=np.sqrt((la.var(ddof=1)+lb.var(ddof=1))/2)
            smd=(la.mean()-lb.mean())/pooled_sd if pooled_sd>0 else np.nan
            matched.append(dict(donor_id=donor.donor_id,disease_status=donor.disease_status,celltype=typ,
                replicate=rep+1,matched_n_per_group=nmatch,mural_percent=ap,endothelial_percent=bp,
                difference_pp=ap-bp,mural_median_UMI=dat.loc[ia,'nCount_RNA'].median(),
                endothelial_median_UMI=dat.loc[ib,'nCount_RNA'].median(),
                log1p_UMI_SMD=smd))
matched=pd.DataFrame(matched);matched.to_csv(OUT/'depth_matched_replicates.csv',index=False)
pd.DataFrame(match_audit).to_csv(OUT/'depth_matching_inclusion.csv',index=False)
md=matched.groupby(['donor_id','disease_status','celltype']).agg(matched_n_per_group=('matched_n_per_group','first'),
    mural_percent=('mural_percent','mean'),endothelial_percent=('endothelial_percent','mean'),
    difference_pp=('difference_pp','mean'),matching_SD_pp=('difference_pp','std'),
    mural_median_UMI=('mural_median_UMI','mean'),endothelial_median_UMI=('endothelial_median_UMI','mean'),
    log1p_UMI_SMD=('log1p_UMI_SMD','mean')).reset_index()
md.to_csv(OUT/'depth_matched_donor_summary.csv',index=False)
ms=[]
for cohort in COHORTS:
    for typ in ['SMC','pericyte']:
        q=md[cohort_mask(md,cohort)&md.celltype.eq(typ)]
        ms.append(dict(cohort=cohort,celltype=typ,donors_n=len(q),mean_difference_pp=q.difference_pp.mean(),
            sd_difference_pp=q.difference_pp.std(ddof=1),min_donor_difference_pp=q.difference_pp.min(),
            max_donor_difference_pp=q.difference_pp.max(),mural_mean_percent=q.mural_percent.mean(),
            endothelial_mean_percent=q.endothelial_percent.mean()))
pd.DataFrame(ms).to_csv(OUT/'depth_matched_cohort_summary.csv',index=False)
audit['finished_utc']=datetime.now(timezone.utc).isoformat()
audit['protocol_unchanged']=sha(protocol)==audit['protocol_sha256']
assert audit['protocol_unchanged']
(OUT/'analysis_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
print(primary[primary.cohort.isin(['all','control'])].to_string(index=False))
print('Depth matching:');print(pd.DataFrame(ms).to_string(index=False))
