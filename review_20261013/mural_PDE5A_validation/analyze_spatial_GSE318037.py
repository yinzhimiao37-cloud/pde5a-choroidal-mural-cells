"""PDE5A阳性壁细胞与内皮的空间条件置换；固定规则见spatial_protocol_locked.md。"""
from pathlib import Path
import hashlib,json,time
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

HERE=Path(__file__).resolve().parent
OUT=HERE/'spatial_GSE318037'
SOURCE=OUT/'joined_all_cells.csv.gz'
N_PERM=9999
SEED=20261013
RADII=np.arange(0,101,10)
NEIGHBOR_RADII=[10,25,50,100]
METRICS=[f'EC_neighbors_{r}um' for r in NEIGHBOR_RADII]+['nearest_EC_um']+[f'nearest_CDF_{r}um' for r in RADII]
PRIMARY=METRICS.index('EC_neighbors_25um')
EC=['artery','vein','choriocapillaris']
started=time.perf_counter()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
protocol_hash=sha(HERE/'spatial_protocol_locked.md')
df=pd.read_csv(SOURCE).rename(columns={'donor_id':'donor','section_id':'section',
    'PDE5A_count':'PDE5A_transcript_count','total_transcripts':'total_transcript_count'})
# 仅使用作者已注释的细胞；保留原始输入以审计作者QC。area_qc还包含核面积完整性，
# 不是本协议纳入标准，不能用它误删合法的细胞边界分割结果。
original_cells_n=len(df)
assert df.author_annotation_available.isin([True,False]).all()
df=df[df.author_annotation_available.eq(True)].copy()
assert df.coordinate_qc.eq(True).all()
assert df.cell_area.gt(0).all()
required=['donor','section','region','x_um','y_um','celltype','PDE5A_transcript_count','total_transcript_count']
assert set(required).issubset(df.columns)
assert df[required].notna().all().all()
assert np.isfinite(df[['x_um','y_um']].to_numpy()).all()
assert df.groupby('section').donor.nunique().eq(1).all()
assert df.section.nunique()==7 and df.donor.nunique()==6
assert df.PDE5A_transcript_count.ge(0).all()
assert np.allclose(df.PDE5A_transcript_count,np.round(df.PDE5A_transcript_count))
all_counts=df.groupby(['donor','section','region','celltype']).size().reset_index(name='cells_n')
all_counts.to_csv(OUT/'all_region_celltype_counts.csv',index=False)

def effect_summary(obs,null):
    mu=float(null.mean());sd=float(null.std(ddof=1))
    if sd<1e-12:return mu,sd,np.nan,np.nan
    nes=float((obs-mu)/sd)
    upper=(1+np.count_nonzero(null>=obs-1e-12))/(len(null)+1)
    lower=(1+np.count_nonzero(null<=obs+1e-12))/(len(null)+1)
    return mu,sd,nes,float(min(1,2*min(upper,lower)))

def bins(x):
    # 依据连续量本身切分三分位，不按PDE5A阳性/阴性结果分箱。
    return pd.qcut(np.asarray(x),q=3,labels=False,duplicates='drop')

def permute_metrics(values,pos,strata,rng):
    npos=int(pos.sum());null=np.zeros((N_PERM,values.shape[1]),dtype=float)
    mobile_cells=0;mobile_positive=0;stratum_rows=[]
    for group in sorted(np.unique(strata)):
        ix=np.flatnonzero(strata==group);k=int(pos[ix].sum());n=len(ix)
        mobile=n>=5 and 0<k<n
        stratum_rows.append({'stratum':str(group),'cells_n':n,'positive_n':k,'exchangeable':mobile})
        if not mobile:
            null += values[ix[pos[ix]]].sum(axis=0)
            continue
        mobile_cells+=n;mobile_positive+=k
        # 分批均匀无放回抽样，避免9999×所有细胞的大内存矩阵。
        for start in range(0,N_PERM,250):
            count=min(250,N_PERM-start)
            sel=np.argpartition(rng.random((count,n)),k-1,axis=1)[:,:k]
            for j in range(values.shape[1]):
                null[start:start+count,j]+=values[ix,j][sel].sum(axis=1)
    return null/npos,mobile_cells,mobile_positive,stratum_rows

sections=[];eligibility=[];strata_audit=[];cell_tables=[]
null_store={};rng=np.random.default_rng(SEED)
for section,region in df[df.region.eq('Choroid')].groupby('section',sort=True):
    donor=region.donor.iloc[0]
    condition='Control' if str(section).startswith('Control') else 'CNV'
    endothelial=region[region.celltype.isin(EC)]
    print('section',section,'choroid cells',len(region),'EC',len(endothelial),flush=True)
    if endothelial.empty:
        eligibility.append(dict(section=section,donor=donor,reason='no_Choroid_endothelium'));continue
    ec_tree=cKDTree(endothelial[['x_um','y_um']].to_numpy())
    region_tree=cKDTree(region[['x_um','y_um']].to_numpy())
    for typ in ['SMC','pericyte']:
        sub=region[region.celltype.eq(typ)].copy()
        if sub.empty:
            eligibility.append(dict(section=section,donor=donor,celltype=typ,reason='no_mural_cells'));continue
        xy=sub[['x_um','y_um']].to_numpy()
        dist=ec_tree.query(xy,k=1)[0]
        neighbors=np.column_stack([ec_tree.query_ball_point(xy,r,return_length=True) for r in NEIGHBOR_RADII])
        density=region_tree.query_ball_point(xy,50,return_length=True)-1
        values=np.column_stack([neighbors,dist,*[(dist<=r).astype(float) for r in RADII]])
        source=sub.copy()
        source['nearest_EC_um']=dist;source['local_Choroid_cells_50um']=density
        for j,r in enumerate(NEIGHBOR_RADII):source[f'EC_neighbors_{r}um']=neighbors[:,j]
        cell_tables.append(source)
        d_bin=np.nan_to_num(bins(density),nan=0).astype(int)
        t_bin=np.nan_to_num(bins(sub.total_transcript_count),nan=0).astype(int)
        for threshold in [1,2]:
            pos=sub.PDE5A_transcript_count.ge(threshold).to_numpy()
            k=int(pos.sum());n=len(sub)
            reason='eligible' if n>=20 and k>=5 and n-k>=5 else 'below_fixed_cell_or_positive_negative_threshold'
            eligibility.append(dict(section=section,donor=donor,condition=condition,celltype=typ,
                PDE5A_min_transcripts=threshold,cells_n=n,positive_n=k,negative_n=n-k,
                endothelial_n=len(endothelial),reason=reason))
            if reason!='eligible':continue
            schemes=['depth_density_stratified','unstratified_within_section'] if threshold==1 else ['depth_density_stratified']
            for scheme in schemes:
                strata=d_bin*3+t_bin if scheme=='depth_density_stratified' else np.zeros(n,dtype=int)
                null,mobile_n,mobile_k,sa=permute_metrics(values,pos,strata,rng)
                key=f'{section}_{typ}_ge{threshold}_{scheme}'
                null_store[key]=null
                observed=values[pos].mean(axis=0)
                negative=values[~pos].mean(axis=0)
                record=dict(section=section,donor=donor,condition=condition,celltype=typ,
                    PDE5A_min_transcripts=threshold,scheme=scheme,cells_n=n,positive_n=k,
                    exchangeable_cells_n=mobile_n,exchangeable_positive_n=mobile_k,
                    exchangeable_positive_fraction=mobile_k/k,null_key=key)
                for j,metric in enumerate(METRICS):
                    record[metric]=float(observed[j]);record[metric+'_negative']=float(negative[j])
                    record[metric+'_null_mean']=float(null[:,j].mean())
                _,nullsd,nes,p=effect_summary(observed[PRIMARY],null[:,PRIMARY])
                record.update(primary_null_sd=nullsd,NES=nes,permutation_p_two_sided=p)
                sections.append(record)
                for st in sa:strata_audit.append({**{a:record[a] for a in ['section','donor','celltype','PDE5A_min_transcripts','scheme']},**st})
pd.DataFrame(eligibility).to_csv(OUT/'spatial_inclusion_audit.csv',index=False)
pd.DataFrame(strata_audit).to_csv(OUT/'permutation_strata_audit.csv',index=False)
pd.concat(cell_tables,ignore_index=True).to_csv(OUT/'Choroid_mural_cell_spatial_source.csv.gz',index=False)
sec=pd.DataFrame(sections)
if sec.empty:
    (OUT/'no_eligible_groups.json').write_text(json.dumps({'reason':'No section and subtype met the frozen eligibility criteria'},indent=2))
    raise SystemExit('No eligible spatial groups; see spatial_inclusion_audit.csv')
sec.to_csv(OUT/'section_spatial_summary.csv',index=False)
np.savez_compressed(OUT/'section_permutation_nulls.npz',**null_store)

# 先在同一donor内合并section，每位供体才是一个汇总点。
donor_rows=[];donor_null={}
for (donor,typ,threshold,scheme),q in sec.groupby(['donor','celltype','PDE5A_min_transcripts','scheme']):
    weights=q.positive_n.to_numpy(dtype=float);weights/=weights.sum()
    obs=np.sum(q[METRICS].to_numpy()*weights[:,None],axis=0)
    null=np.sum(np.stack([null_store[k] for k in q.null_key])*weights[:,None,None],axis=0)
    key=f'{donor}_{typ}_ge{threshold}_{scheme}'
    donor_null[key]=null
    mu,sd,nes,p=effect_summary(obs[PRIMARY],null[:,PRIMARY])
    r=dict(donor=donor,condition=q.condition.iloc[0],celltype=typ,PDE5A_min_transcripts=threshold,scheme=scheme,
        sections_n=len(q),section_ids=';'.join(q.section),positive_n=int(q.positive_n.sum()),
        cells_n=int(q.cells_n.sum()),exchangeable_positive_fraction=q.exchangeable_positive_n.sum()/q.positive_n.sum(),
        null_key=key,NES=nes,permutation_p_two_sided=p,primary_null_sd=sd)
    for j,metric in enumerate(METRICS):
        r[metric]=obs[j];r[metric+'_null_mean']=null[:,j].mean()
    donor_rows.append(r)
don=pd.DataFrame(donor_rows);don.to_csv(OUT/'donor_spatial_summary.csv',index=False)
np.savez_compressed(OUT/'donor_permutation_nulls.npz',**donor_null)
cohort_rows=[];curve_rows=[]
for (typ,threshold,scheme),base in don.groupby(['celltype','PDE5A_min_transcripts','scheme']):
    for cohort in ['all','Control','CNV','exclude_possible_shared_Control3']:
        q=base if cohort=='all' else (base[~base.donor.eq('Control3')] if cohort=='exclude_possible_shared_Control3' else base[base.condition.eq(cohort)])
        if q.empty:continue
        obs=q[METRICS].mean().to_numpy()
        null=np.stack([donor_null[k] for k in q.null_key]).mean(axis=0)
        mu,sd,nes,p=effect_summary(obs[PRIMARY],null[:,PRIMARY])
        cohort_rows.append(dict(cohort=cohort,celltype=typ,PDE5A_min_transcripts=threshold,scheme=scheme,
            donors_n=len(q),donor_ids=';'.join(map(str,q.donor)),sections_n=int(q.sections_n.sum()),
            positive_cells_n=int(q.positive_n.sum()),mean_neighbors_25um=obs[PRIMARY],
            sd_donor_neighbors_25um=q['EC_neighbors_25um'].std(ddof=1),null_mean=mu,null_sd=sd,NES=nes,
            permutation_p_two_sided=p if cohort=='all' else np.nan,
            mean_nearest_distance_um=obs[4],sd_donor_nearest_distance_um=q['nearest_EC_um'].std(ddof=1),
            exchangeable_positive_fraction=q.exchangeable_positive_fraction.mean()))
        for radius in RADII:
            metric=f'nearest_CDF_{radius}um';j=METRICS.index(metric)
            curve_rows.append(dict(cohort=cohort,celltype=typ,PDE5A_min_transcripts=threshold,scheme=scheme,
                radius_um=radius,donors_n=len(q),observed_fraction=obs[j],sd_across_donors=q[metric].std(ddof=1),
                null_mean_fraction=null[:,j].mean(),null_lower_025=np.quantile(null[:,j],.025),
                null_upper_975=np.quantile(null[:,j],.975)))
co=pd.DataFrame(cohort_rows);co['BH_q_primary_family']=np.nan
for (threshold,scheme),q in co[co.cohort.eq('all')].groupby(['PDE5A_min_transcripts','scheme']):
    p=q.permutation_p_two_sided.dropna().sort_values()
    if len(p):
        # 冻结家族为两个亚型；不可评价亚型仍占一个预定检验位置（相当于p=1）。
        adjusted=np.minimum.accumulate((p.to_numpy()*2/np.arange(1,len(p)+1))[::-1])[::-1]
        co.loc[p.index,'BH_q_primary_family']=np.minimum(adjusted,1)
co.to_csv(OUT/'cohort_spatial_summary.csv',index=False)
pd.DataFrame(curve_rows).to_csv(OUT/'distance_gradient_source.csv',index=False)
audit={'input_sha256':sha(SOURCE),'protocol_sha256':protocol_hash,
    'input_segmented_cells_n':original_cells_n,'author_annotated_cells_n':len(df),
    'protocol_unchanged':protocol_hash==sha(HERE/'spatial_protocol_locked.md'),
    'permutations':N_PERM,'seed':SEED,'main_region':'Choroid','unit':'donor',
    'metric':'Mean endothelial centroid neighbours within 25 micrometres per PDE5A-positive mural cell',
    'coordinates':'cell centroids; not membrane contact','elapsed_seconds':time.perf_counter()-started,
    'multiple_testing':'Two planned mural subtypes in each positivity-threshold and scheme family, including unevaluable tests in family size; only all-donor P values inferential',
    'limitations':'Control/CNV, not CSCR/tadalafil; conditional label exchangeability, residual local architecture and technical confounding possible'}
(OUT/'spatial_analysis_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
print(co[co.cohort.eq('all')].to_string(index=False))
