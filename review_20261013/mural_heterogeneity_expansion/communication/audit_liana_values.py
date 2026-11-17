"""从共享原始RNA矩阵独立复算100个LIANA评分，核查分组/复合体/ID一致性。"""
from pathlib import Path
import json
import anndata as ad
import pandas as pd
import numpy as np

HERE = Path(__file__).resolve().parent
root = HERE.parents[2]
a = ad.read_h5ad(root/'review_20261013/mural_heterogeneity_extension/scRNA/GSE183320_mural_EC_macrophage_normalized.h5ad')
obs = a.obs.copy()
obs['donor_id'] = obs.donor_id.astype(int)
obs['PDE5A'] = a[:,'PDE5A'].X.toarray().ravel()
s = pd.read_csv(HERE/'results/LIANA_all_donor_scores.csv.gz')
s = s.loc[s.variant.eq('high_vs_low_positive__unmatched') & s.lrs_to_keep]
pick = s.sample(min(100,len(s)),random_state=20261013)
rows=[]
for _, r in pick.iterrows():
    dm = obs.loc[obs.donor_id.eq(r.donor_id)]
    mt = dm.loc[dm.celltype.eq(r.mural_type) & dm.PDE5A.gt(0)]
    median = mt.PDE5A.median()
    sm = mt.loc[mt.PDE5A.gt(median) if r.state=='high' else mt.PDE5A.le(median)]
    types = ['artery','vein','choriocapillaris','ec_unknown'] if r.target=='endothelial' else ['macrophage-inf','macrophage-res']
    rm = dm.loc[dm.celltype.isin(types)]
    ligand = r.ligand_complex.split('_')
    receptor = r.receptor_complex.split('_')
    # LIANA geometric_mean按复合体最小亚单位均值组装，再对配体/受体均值取几何均值。
    lm = float(np.asarray(a[sm.index,ligand].X.mean(axis=0)).ravel().min())
    rmv = float(np.asarray(a[rm.index,receptor].X.mean(axis=0)).ravel().min())
    expected = np.sqrt(lm*rmv)
    rows.append(dict(donor_id=int(r.donor_id),mural_type=r.mural_type,state=r.state,
        target=r.target,ligand=r.ligand_complex,receptor=r.receptor_complex,
        liana_score=r.lr_gmeans,independent_score=expected,absolute_error=abs(expected-r.lr_gmeans)))
q=pd.DataFrame(rows)
q.to_csv(HERE/'results/independent_score_audit.csv',index=False)
assert q.absolute_error.max()<3e-6

stats=pd.read_csv(HERE/'results/donor_LR_statistics_all_families.csv')
paired=pd.read_csv(HERE/'results/donor_paired_LR_differences.csv.gz')
assert not paired.duplicated(['variant','donor_id','mural_type','target','ligand_complex','receptor_complex']).any()
assert (stats.donors_n>=5).all()
primary=stats[(stats.variant=='high_vs_low_positive__unmatched')&(stats.cohort=='all')]
audit={'score_rows_independently_checked':len(q),'max_score_absolute_error':float(q.absolute_error.max()),
       'no_duplicate_donor_LR_pairs':True,'all_inferential_tests_donors_ge5':True,
       'primary_tested_LR':len(primary),'primary_BH_significant':int((primary.q_BH<.05).sum()),
       'primary_positive_candidates':int(primary.candidate_support.sum())}
(HERE/'results/QA_independent_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
print(json.dumps(audit,indent=2))

# 完整列出12个固定家族，包括不可评价的家族。
family=stats.groupby(['variant','cohort']).agg(tested_LR=('p_value','size'),candidate_LR=('candidate_support','sum'),
    min_donors=('donors_n','min'),max_donors=('donors_n','max'),min_q=('q_BH','min'))
idx=pd.MultiIndex.from_product([sorted(paired.variant.unique()),['all','control_only','exclude_known_shared']],names=['variant','cohort'])
family=family.reindex(idx)
family[['tested_LR','candidate_LR']]=family[['tested_LR','candidate_LR']].fillna(0).astype(int)
family['status']=np.where(family.tested_LR>0,'evaluated','no LR with >=5 eligible donors')
family.reset_index().to_csv(HERE/'results/multiple_testing_family_summary.csv',index=False)
