"""一次分块读取 GSE183320；为壁细胞分析和通讯分析共享全基因稀疏矩阵。

运行：py review_20261013/mural_heterogeneity_extension/scRNA/prepare_expression.py
依赖：pandas numpy scipy anndata
已标准化数据不再 log、不转成伪造的整数 UMI。
"""
from pathlib import Path
import csv, gzip, json, hashlib
import numpy as np
import pandas as pd
from scipy import sparse
import anndata as ad

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[2]
META = PROJECT/'data/raw/GSE183320_all_metadata.csv.gz'
EXPR = PROJECT/'data/raw/GSE183320_all_normalized_expression_data.csv.gz'
OUT = HERE/'GSE183320_mural_EC_macrophage_normalized.h5ad'
KEEP = ['SMC','pericyte','artery','vein','choriocapillaris','ec_unknown','macrophage-inf','macrophage-res']

def main():
    if OUT.exists():
        print(f'Cache exists: {OUT}', flush=True); return
    meta=pd.read_csv(META)
    meta['barcode']=meta['barcode'].astype(str).str.strip()
    meta=meta.loc[meta.celltype.isin(KEEP)].copy().set_index('barcode')
    meta['donor_id']=meta.library.str.extract(r'donor_(\d+)',expand=False).astype(str)
    with gzip.open(EXPR,'rt') as f: header=next(csv.reader(f))
    # R 导出把 -1_ 改成 .1_；只作精确格式复原，之后必须一对一检查。
    fixed=[s.strip().replace('.1_','-1_') for s in header[1:]]
    assert len(set(fixed))==len(fixed)
    mapping=dict(zip(fixed,header[1:]))
    missing=set(meta.index)-set(mapping)
    assert not missing, missing
    cols=[mapping[b] for b in meta.index]
    blocks=[]; names=[]
    for i,df in enumerate(pd.read_csv(EXPR,usecols=[header[0]]+cols,chunksize=128)):
        gene=df[header[0]].astype(str).str.strip()
        x=df[cols].to_numpy(dtype=np.float32)
        if not np.isfinite(x).all() or np.min(x)<0:
            raise ValueError('Source expression must be finite/nonnegative')
        blocks.append(sparse.csr_matrix(x)); names.extend(gene)
        if i%20==0: print(f'Read {len(names)} genes',flush=True)
    assert len(set(names))==len(names), 'Duplicated gene names require explicit resolution'
    X=sparse.vstack(blocks,format='csr').T.tocsr()
    obj=ad.AnnData(X,obs=meta,var=pd.DataFrame(index=names))
    obj.uns['expression_provenance']='Author normalized matrix; no re-normalization or second log; not raw counts'
    # 所有字符串列统一，避免 HDF5 object 混合类型写入失败。
    for c in obj.obs:
        if obj.obs[c].dtype==object: obj.obs[c]=obj.obs[c].fillna('NA').astype(str)
    obj.write_h5ad(OUT,compression='gzip')
    audit={'source_expression':str(EXPR),'source_metadata':str(META),'n_cells':obj.n_obs,'n_genes':obj.n_vars,
           'celltypes':obj.obs.celltype.value_counts().to_dict(),'expression_transform':'none','matrix_is_raw_counts':False}
    (HERE/'expression_cache_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    print(json.dumps(audit,indent=2),flush=True)

if __name__=='__main__': main()
