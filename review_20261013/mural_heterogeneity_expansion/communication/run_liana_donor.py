#!/usr/bin/env python
"""GSE183320 壁细胞 PDE5A 分层 LIANA 分析；独立统计单位为供体。

运行：py run_liana_donor.py --stage all
本地 .runtime 只为避免改动原项目环境；一般用户直接 pip install -r requirements.txt。
不把零/缺失评分填补为支持通讯，也不把细胞数当作供体数。
"""
from __future__ import annotations
import argparse, hashlib, itertools, json, sys, warnings
from pathlib import Path

HERE = Path(__file__).resolve().parent
if (HERE / '.runtime').exists():
    sys.path.insert(0, str(HERE / '.runtime'))
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from statsmodels.stats.multitest import multipletests
import liana as li

ROOT = HERE.parents[2]
DEFAULT_INPUT = ROOT / 'review_20261013/mural_heterogeneity_extension/scRNA/GSE183320_mural_EC_macrophage_normalized.h5ad'
OUT = HERE / 'results'
SEED = 20261013
SENDER_TYPES = ('SMC', 'pericyte')
RECEIVER_TYPES = ('endothelial', 'macrophage')
KEYS = ['mural_type', 'target', 'ligand_complex', 'receptor_complex']


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def smd(high, low, column):
    """在 log1p 技术指标上计算组间标准化均值差，不输出未经校正的生物学P。"""
    a, b = np.log1p(high[column].to_numpy(float)), np.log1p(low[column].to_numpy(float))
    if min(len(a), len(b)) < 2:
        return np.nan
    sd = np.sqrt((np.var(a, ddof=1) + np.var(b, ddof=1)) / 2)
    return (a.mean() - b.mean()) / sd if sd > 0 else np.nan


def depth_match(high, low, rng):
    """事前固定：pooled nFeature 四分位分箱，箱内等量、无放回抽样一次。"""
    both = pd.concat([high.assign(_arm='high'), low.assign(_arm='low')])
    if both.empty:
        return high, low
    both['_bin'] = pd.qcut(both.nFeature_RNA, q=4, duplicates='drop', labels=False).fillna(0)
    hi, lo = [], []
    for _, part in both.groupby('_bin', observed=True):
        a = part.index[part._arm == 'high'].to_numpy()
        b = part.index[part._arm == 'low'].to_numpy()
        k = min(len(a), len(b))
        if k:
            hi.extend(rng.choice(a, k, replace=False))
            lo.extend(rng.choice(b, k, replace=False))
    return high.loc[hi], low.loc[lo]


def infer(input_path):
    OUT.mkdir(parents=True, exist_ok=True)
    a = ad.read_h5ad(input_path)
    if not sparse.issparse(a.X):
        a.X = sparse.csr_matrix(a.X)
    a.X = a.X.astype(np.float32).tocsr()
    if a.X.data.min(initial=0) < 0:
        raise ValueError('检测到负值：不能使用 integrated/scale.data 进行此分析。')
    if not a.var_names.is_unique or not a.obs_names.is_unique:
        raise ValueError('基因/细胞ID不唯一，需审计后解决，不能静默去重。')
    if 'PDE5A' not in a.var_names:
        raise ValueError('输入缺少 PDE5A。')
    obs = a.obs.copy()
    obs['donor_id'] = obs.donor_id.astype(int)
    obs['PDE5A_value'] = a[:, 'PDE5A'].X.toarray().ravel()
    obs['receiver'] = obs.celltype.map({
        'artery': 'endothelial', 'vein': 'endothelial',
        'choriocapillaris': 'endothelial', 'ec_unknown': 'endothelial',
        'macrophage-inf': 'macrophage', 'macrophage-res': 'macrophage',
        'macrophage': 'macrophage'})
    resource = li.resource.select_resource('consensus')
    # PDE5A 不是分泌配体；防止任何资源未来版本将其作为候选输入。
    bad = resource.apply(lambda r: 'PDE5A' in (str(r.ligand)+'_'+str(r.receptor)).split('_'), axis=1)
    resource = resource.loc[~bad].copy()
    resource.to_csv(OUT / 'LIANA_consensus_resource.csv', index=False)
    gene_set = set(a.var_names)
    resource['all_subunits_available'] = resource.apply(
        lambda r: set((r.ligand+'_'+r.receptor).split('_')).issubset(gene_set), axis=1)
    resource.to_csv(OUT / 'resource_coverage_audit.csv', index=False)
    resource = resource.loc[resource.all_subunits_available, ['ligand', 'receptor']]
    a = a[:, a.var_names != 'PDE5A'].copy()
    # 资源数据、PDE5A分组与技术指标分别保留；原RNA矩阵不再次归一化。
    audits, depths, definitions, all_res = [], [], [], []
    for donor in sorted(obs.donor_id.unique()):
        dm = obs.loc[obs.donor_id == donor].copy()
        disease = str(dm.disease_status.iloc[0])
        donor_groups = {}
        for typ in SENDER_TYPES:
            mt = dm.loc[dm.celltype == typ]
            pos = mt.loc[mt.PDE5A_value > 0]
            med = float(pos.PDE5A_value.median()) if len(pos) else np.nan
            donor_groups[typ] = {
                'high': pos.loc[pos.PDE5A_value > med],
                'low_positive': pos.loc[pos.PDE5A_value <= med],
                'undetected': mt.loc[mt.PDE5A_value == 0]}
            for state, cells in donor_groups[typ].items():
                definitions.append(dict(donor_id=donor, disease_status=disease, mural_type=typ,
                    state=state, positive_median=med, cells_n=len(cells),
                    median_PDE5A=float(cells.PDE5A_value.median()) if len(cells) else np.nan))
        for comp in ('low_positive', 'undetected'):
            for matching in ('unmatched', 'nFeature_matched'):
                variant = f'high_vs_{comp}__{matching}'
                labels = pd.Series('', index=dm.index, dtype=object)
                pairs = []
                for rec in RECEIVER_TYPES:
                    ridx = dm.index[dm.receiver == rec]
                    if len(ridx) >= 20:
                        labels.loc[ridx] = rec
                for ti, typ in enumerate(SENDER_TYPES):
                    high, low = donor_groups[typ]['high'], donor_groups[typ][comp]
                    if matching == 'nFeature_matched':
                        # 唯一固定种子与donor/type/comparison键，不根据结果尝试新种子。
                        rng = np.random.default_rng(SEED + donor*100 + ti*10 + int(comp=='undetected'))
                        high, low = depth_match(high, low, rng)
                    depths.append(dict(donor_id=donor, disease_status=disease, mural_type=typ,
                        variant=variant, high_n=len(high), comparator_n=len(low),
                        nCount_log_SMD=smd(high, low, 'nCount_RNA'),
                        nFeature_log_SMD=smd(high, low, 'nFeature_RNA'),
                        high_median_nCount=high.nCount_RNA.median(),
                        comparator_median_nCount=low.nCount_RNA.median(),
                        high_median_nFeature=high.nFeature_RNA.median(),
                        comparator_median_nFeature=low.nFeature_RNA.median()))
                    enough = min(len(high), len(low)) >= 10
                    for rec in RECEIVER_TYPES:
                        nr = int((dm.receiver == rec).sum())
                        valid = enough and nr >= 20
                        audits.append(dict(donor_id=donor, disease_status=disease, mural_type=typ,
                            target=rec, variant=variant, high_n=len(high), comparator_n=len(low),
                            receiver_n=nr, eligible=valid,
                            reason='eligible' if valid else 'sender<10 or receiver<20'))
                        if valid:
                            pairs.extend([{'source': typ+'__'+state, 'target': rec}
                                          for state in ('high', comp)])
                    if enough:
                        labels.loc[high.index] = typ+'__high'
                        labels.loc[low.index] = typ+'__'+comp
                if not pairs:
                    continue
                kept = labels.index[labels != '']
                sub = a[kept].copy()
                sub.obs['lr_group'] = pd.Categorical(labels.loc[kept])
                # n_perms=None明确不计算细胞层面的显著性。
                res = li.mt.geometric_mean(sub, groupby='lr_group', resource=resource,
                    expr_prop=0.10, min_cells=10, groupby_pairs=pd.DataFrame(pairs),
                    n_perms=None, return_all_lrs=False, use_raw=False, inplace=False,
                    verbose=False, seed=SEED)
                # LIANA 1.10.0 的 return_all_lrs=True 分支在复合体汇总前去重，
                # 独立复算证实会错取首个亚单位。因此用False正确组装复合体，
                # 再与完整资源左连接。未通过支持门槛/未表达者=NaN，不能赋零。
                res['lrs_to_keep'] = True
                full = pd.DataFrame(pairs).merge(resource.rename(columns={
                    'ligand':'ligand_complex','receptor':'receptor_complex'}),how='cross')
                res = full.merge(res,on=['source','target','ligand_complex','receptor_complex'],
                    how='left',validate='one_to_one')
                res['lrs_to_keep'] = res.lrs_to_keep.fillna(False).astype(bool)
                res['donor_id'], res['disease_status'], res['variant'] = donor, disease, variant
                res['mural_type'] = res.source.str.split('__').str[0]
                res['state'] = res.source.str.split('__').str[1]
                all_res.append(res)
                print(f'donor={donor}, {variant}, {len(res)} LR rows', flush=True)
    pd.DataFrame(definitions).to_csv(OUT / 'PDE5A_group_definition_audit.csv', index=False)
    pd.DataFrame(audits).to_csv(OUT / 'donor_group_inclusion.csv', index=False)
    pd.DataFrame(depths).to_csv(OUT / 'sender_depth_balance.csv', index=False)
    if not all_res:
        raise RuntimeError('固定细胞数规则下没有可评价sender–receiver，未输出假零结果。')
    scores = pd.concat(all_res, ignore_index=True)
    scores.to_csv(OUT / 'LIANA_all_donor_scores.csv.gz', index=False)
    audit = dict(dataset='GSE183320',input=str(input_path),input_sha256=sha256(input_path),
        cells_n=int(a.n_obs), donor_n=int(obs.donor_id.nunique()),
        genes_n=int(a.n_vars), liana=li.__version__, pandas=pd.__version__,
        numpy=np.__version__, anndata=ad.__version__, random_seed=SEED,
        resource_interactions_available=int(len(resource)),
        resource_sha256=sha256(OUT/'LIANA_consensus_resource.csv'),
        score_rows=int(len(scores)), cell_level_permutation_pvalues_used=False)
    (OUT/'analysis_audit.json').write_text(json.dumps(audit, indent=2), encoding='utf-8')


def donor_tests(deltas, seed=SEED):
    """按实际供体集合分批。符号翻转是供体差值对称性假设下的推断，不是随机试验。"""
    rows = []
    for key, part in deltas.groupby(KEYS, observed=True, sort=True):
        d = part.delta.to_numpy(float)
        n = len(d)
        # 每位供体一个且仅一个配对值。
        assert part.donor_id.nunique() == n
        if n < 5:
            continue
        rng = np.random.default_rng(seed)
        effect = float(d.mean())
        if n <= 15:
            bits = ((np.arange(2**n, dtype=np.uint32)[:,None] >> np.arange(n)) & 1)
            signs = bits.astype(float)*2-1
            p = float((np.abs(signs @ d / n) >= abs(effect)-1e-12).mean())
            how = 'exact donor sign-flip'
        else:
            # 分块避免大型中间矩阵；不是把99999个样本当成供体。
            exceed, done = 0, 0
            while done < 99999:
                k = min(5000, 99999-done)
                signs = rng.choice(np.array([-1., 1.]), size=(k,n))
                exceed += int((np.abs(signs @ d/n) >= abs(effect)-1e-12).sum())
                done += k
            p = (exceed+1)/100000
            how = '99999 donor sign-flips plus-one'
        boot = np.mean(d[rng.integers(0,n,size=(9999,n))],axis=1)
        ci = np.quantile(boot,[0.025,0.975])
        loo = (d.sum()-d)/(n-1)
        row = dict(zip(KEYS,key))
        row.update(donors_n=n, donor_ids=';'.join(map(str, sorted(part.donor_id.astype(int)))),
            mean_delta=effect, sd_delta=float(d.std(ddof=1)),
            bootstrap95_low=float(ci[0]), bootstrap95_high=float(ci[1]),
            positive_fraction=float(np.mean(d>0)), p_value=p, test=how,
            loo_mean_min=float(loo.min()), loo_mean_max=float(loo.max()),
            loo_positive_all=bool(np.all(loo>0)),
            mean_high=float(part.high.mean()), mean_comparator=float(part.comparator.mean()))
        rows.append(row)
    return pd.DataFrame(rows)


def summarize():
    scores = pd.read_csv(OUT/'LIANA_all_donor_scores.csv.gz')
    # 仅纳入LIANA通过门槛者；完整资源中补回的缺失结果绝不能当成真实零。
    supported = scores.loc[scores.lrs_to_keep.astype(str).str.lower().eq('true')].copy()
    keepkeys = ['variant','donor_id','disease_status']+KEYS
    main = supported.loc[supported.state == 'high', keepkeys+['lr_gmeans']].rename(columns={'lr_gmeans':'high'})
    comp = supported.loc[supported.state != 'high', keepkeys+['lr_gmeans']].rename(columns={'lr_gmeans':'comparator'})
    paired = main.merge(comp, on=keepkeys, how='inner', validate='one_to_one')
    paired['delta'] = paired.high-paired.comparator
    paired.to_csv(OUT/'donor_paired_LR_differences.csv.gz',index=False)
    allsummary = []
    for variant in sorted(paired.variant.unique()):
        vp = paired.loc[paired.variant == variant]
        for cohort in ('all','control_only','exclude_known_shared'):
            sub = vp
            if cohort == 'control_only':
                sub = sub.loc[sub.disease_status.eq('control')]
            elif cohort == 'exclude_known_shared':
                sub = sub.loc[~sub.donor_id.isin([1,2,11,12])]
            summary = donor_tests(sub)
            if summary.empty:
                continue
            summary['q_BH'] = multipletests(summary.p_value,method='fdr_bh')[1]
            summary['variant'], summary['cohort'] = variant, cohort
            summary['candidate_support'] = (summary.q_BH<.05)&(summary.mean_delta>0)&(summary.positive_fraction>=.7)
            allsummary.append(summary)
            print(f'{variant}, {cohort}: {len(summary)} tests; {summary.candidate_support.sum()} candidate rows',flush=True)
    if allsummary:
        summary = pd.concat(allsummary,ignore_index=True)
        summary.to_csv(OUT/'donor_LR_statistics_all_families.csv',index=False)
        primary = summary.loc[(summary.variant=='high_vs_low_positive__unmatched')&(summary.cohort=='all')]
        primary.sort_values(['q_BH','mean_delta'],ascending=[True,False]).to_csv(OUT/'primary_all_tested_LR.csv',index=False)
        primary.loc[primary.candidate_support].sort_values('q_BH').to_csv(OUT/'primary_candidate_LR.csv',index=False)
        family = summary.groupby(['variant','cohort'],observed=True).agg(
            tested_LR=('p_value','size'), candidate_LR=('candidate_support','sum'),
            min_donors=('donors_n','min'), max_donors=('donors_n','max'),
            min_q=('q_BH','min'))
        # 明确保留因有效供体不足而没有检验的家族，不能把它们静默删除。
        all_families = pd.MultiIndex.from_product([
            sorted(paired.variant.unique()), ['all','control_only','exclude_known_shared']],
            names=['variant','cohort'])
        family = family.reindex(all_families)
        family[['tested_LR','candidate_LR']] = family[['tested_LR','candidate_LR']].fillna(0).astype(int)
        family['status'] = np.where(family.tested_LR>0,'evaluated','no LR with >=5 eligible donors')
        family = family.reset_index()
        family.to_csv(OUT/'multiple_testing_family_summary.csv',index=False)
    else:
        (OUT/'NO_EVALUABLE_DONOR_TESTS.txt').write_text('No LR contrast had >=5 eligible donors; do not substitute cell-level tests.',encoding='utf-8')


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=DEFAULT_INPUT)
    parser.add_argument('--stage',choices=['infer','summarize','all'],default='all')
    args=parser.parse_args()
    warnings.filterwarnings('ignore',category=FutureWarning)
    if args.stage in ('infer','all'):
        infer(args.input)
    if args.stage in ('summarize','all'):
        summarize()
