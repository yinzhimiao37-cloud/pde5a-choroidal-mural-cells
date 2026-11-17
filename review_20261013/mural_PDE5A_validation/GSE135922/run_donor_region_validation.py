"""GSE135922: separate experimental cluster namespaces, marker-first annotation.

Usage: python run_donor_region_validation.py --stage markers
       python run_donor_region_validation.py --stage quantify
The quantification stage requires a saved, reviewed annotation.json.
Original inputs and previous outputs are never edited.
"""
import argparse, csv, gzip, hashlib, json, platform, re, shlex, time
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PROTOCOL = json.loads((HERE / 'fixed_protocol.json').read_text(encoding='utf-8'))
RAW = ROOT / 'data/raw'
FILES = sorted(RAW.glob('GSM4037*expression.tsv.gz'))
assert len(FILES) == 14
MARKERS = PROTOCOL['marker_panel']
MURAL = ['RGS5', 'ACTA2', 'PDGFRB', 'NOTCH3']

def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda:f.read(2**20), b''): h.update(b)
    return h.hexdigest()

def info(p):
    m = re.fullmatch(r'(GSM\d+)_(macula|peripheral)_donor_(\d+)(_enriched)?_expression.tsv.gz', p.name)
    assert m, p
    gsm, region, donor, enriched = m.groups()
    return {'gsm':gsm, 'donor':int(donor), 'region':region,
            'experiment':'CD31_enriched' if enriched else 'unenriched',
            'health':PROTOCOL['disease_map'][donor], 'source_file':p.name}

def read_selected(p, genes):
    # 文件扩展名虽为tsv，实际是R write.table空格分隔，含额外未命名行号。
    # 显式补足行号列，防止usecols导致列错位；只读取必要列。
    with gzip.open(p,'rt') as f:
        header = shlex.split(f.readline())
    wanted = ['barcode','library','final_cluster_labels'] + genes
    missing = sorted(set(wanted)-set(header))
    if missing: raise ValueError(f'{p.name}: missing {missing}')
    frames = []
    for x in pd.read_csv(p, sep=r'\s+', skiprows=1, header=None,
                         names=['row_number']+header, usecols=wanted,
                         chunksize=500, low_memory=False):
        if x[genes].isna().any().any(): raise ValueError('Missing expression')
        if (x[genes] < 0).any().any(): raise ValueError('Unexpected negative normalized values')
        x['final_cluster_labels'] = x['final_cluster_labels'].astype(int)
        frames.append(x)
    x = pd.concat(frames, ignore_index=True)
    assert not x.barcode.duplicated().any(), 'duplicate barcode within region'
    assert x.library.nunique() == 1
    return x

def memory_mb():
    try:
        import psutil
        mem = psutil.Process().memory_info()
        return {'rss_MiB':mem.rss/2**20, 'peak_wset_MiB':getattr(mem,'peak_wset',0)/2**20}
    except ImportError:
        # Windows自带API提供进程峰值工作集，不额外安装包。
        try:
            import ctypes
            from ctypes import wintypes
            class PMC(ctypes.Structure):
                _fields_=[('cb',wintypes.DWORD),('PageFaultCount',wintypes.DWORD)]+[(k,ctypes.c_size_t) for k in ['PeakWorkingSetSize','WorkingSetSize','QuotaPeakPagedPoolUsage','QuotaPagedPoolUsage','QuotaPeakNonPagedPoolUsage','QuotaNonPagedPoolUsage','PagefileUsage','PeakPagefileUsage']]
            pm=PMC();pm.cb=ctypes.sizeof(pm)
            ctypes.windll.kernel32.GetCurrentProcess.restype=wintypes.HANDLE
            ctypes.windll.psapi.GetProcessMemoryInfo.argtypes=[wintypes.HANDLE,ctypes.POINTER(PMC),wintypes.DWORD]
            ctypes.windll.psapi.GetProcessMemoryInfo.restype=wintypes.BOOL
            ok=ctypes.windll.psapi.GetProcessMemoryInfo(ctypes.windll.kernel32.GetCurrentProcess(),ctypes.byref(pm),pm.cb)
            if ok:return {'rss_MiB':pm.WorkingSetSize/2**20,'peak_wset_MiB':pm.PeakWorkingSetSize/2**20}
        except Exception:pass
        return {'unavailable':'process memory instrumentation failed'}

def markers():
    rows=[]
    for p in FILES:
        meta=info(p); x=read_selected(p,MARKERS)
        for cl,g in x.groupby('final_cluster_labels'):
            rows.append(meta | {'cluster':int(cl),'cells_n':len(g)} |
                        {f'{v}_sum':float(g[v].sum()) for v in MARKERS})
        print('markers',p.name,len(x),flush=True)
    sample=pd.DataFrame(rows); sample.to_csv(HERE/'marker_sample_cluster_sums.csv',index=False)
    q=sample.groupby(['experiment','cluster'],as_index=False)[['cells_n']+[f'{v}_sum' for v in MARKERS]].sum()
    for v in MARKERS: q[v]=q[f'{v}_sum']/q.cells_n
    q=q.drop(columns=[f'{v}_sum' for v in MARKERS])
    for exp,ix in q.groupby('experiment').groups.items():
        z=q.loc[ix,MURAL]; scaled=z/z.max(axis=0)
        q.loc[ix,'mural_marker_score']=scaled.mean(axis=1)
        for v in ['RGS5','ACTA2']:q.loc[ix,f'{v}_rank']=z[v].rank(ascending=False,method='min')
    q.to_csv(HERE/'marker_cluster_means_by_experiment.csv',index=False)
    top=q.loc[q.groupby('experiment').mural_marker_score.idxmax()]
    top.to_csv(HERE/'provisional_mural_annotation_marker_only.csv',index=False)
    print(top[['experiment','cluster','cells_n','mural_marker_score','RGS5_rank','ACTA2_rank']].to_string(index=False))

def summarize(d, unit, region='pooled'):
    rows=[]
    for threshold in [20,1,5,10,30,50]:
      for label in ['all','no_noted_ophthalmologic_disease','unenriched','CD31_enriched']:
        z=d[d.cells_n.ge(threshold)].copy()
        if label=='no_noted_ophthalmologic_disease':z=z[z.health.eq(label)]
        elif label in ['unenriched','CD31_enriched']:z=z[z.experiment.eq(label)]
        for anno in ['mural-like']:
          v=z[z.annotation.eq(anno)]; a=v.positivity_percent
          rows.append({'unit':unit,'region':region,'minimum_cells':threshold,'restriction':label,'annotation':anno,
             'donors_n':int(v.donor.nunique()),'observations_n':len(v),'total_cells':int(v.cells_n.sum()),
             'total_positive_cells':int(v.positive_n.sum()),'mean_percent':a.mean(),'sample_sd_percent':a.std(ddof=1),
             'min_percent':a.min(),'max_percent':a.max(),'donors_with_positive_cells':int(v.loc[v.positive_n.gt(0),'donor'].nunique()),
             'pooled_cell_percent':100*v.positive_n.sum()/v.cells_n.sum() if len(v) else np.nan,
             'donor_ids':';'.join(map(str,sorted(v.donor.unique())))})
    return pd.DataFrame(rows)

def quantify():
    annotation=json.loads((HERE/'annotation.json').read_text(encoding='utf-8'))
    rows=[]; counts=[]; marker_checks=[]
    for p in FILES:
        meta=info(p); x=read_selected(p,['PDE5A']+MARKERS)
        cl=annotation['mural_clusters'][meta['experiment']]
        g=x[x.final_cluster_labels.eq(cl)]
        # 零细胞组合保留在审计中，但不作为0%供体点。
        rows.append(meta|{'annotation':'mural-like','original_cluster':cl,'cells_n':len(g),
            'positive_n':int(g.PDE5A.gt(0).sum()),'expression_sum':float(g.PDE5A.sum())})
        counts.append(meta|{'all_cells_n':len(x),'clusters_present':';'.join(map(str,sorted(x.final_cluster_labels.unique())))})
        for c,gg in x.groupby('final_cluster_labels'):
            marker_checks.append(meta|{'cluster':int(c),'cells_n':len(gg)}|{v:float(gg[v].mean()) for v in MARKERS+['PDE5A']})
        print('quantify',p.name,'mural_cells',len(g),flush=True)
    region=pd.DataFrame(rows)
    region['positivity_percent']=np.where(region.cells_n.gt(0),100*region.positive_n/region.cells_n,np.nan)
    region['mean_published_normalized_expression']=region.expression_sum/region.cells_n.replace(0,np.nan)
    region['eligible_primary_20cells']=region.cells_n.ge(20)
    region.to_csv(HERE/'PDE5A_mural_donor_region.csv',index=False)
    keys=['donor','experiment','health','annotation']
    donor=region.groupby(keys,as_index=False)[['cells_n','positive_n','expression_sum']].sum()
    donor['positivity_percent']=100*donor.positive_n/donor.cells_n.replace(0,np.nan)
    donor['mean_published_normalized_expression']=donor.expression_sum/donor.cells_n.replace(0,np.nan)
    donor['eligible_primary_20cells']=donor.cells_n.ge(20)
    donor.to_csv(HERE/'PDE5A_mural_donor_pooled_regions.csv',index=False)
    ds=summarize(donor,'donor_pooled_regions'); ds.to_csv(HERE/'PDE5A_mural_donor_summary_thresholds.csv',index=False)
    rs=pd.concat([summarize(z,'donor_region',reg) for reg,z in region.groupby('region')],ignore_index=True)
    rs.to_csv(HERE/'PDE5A_mural_region_summary_thresholds.csv',index=False)
    lod=[]
    for restriction in ['all','no_noted_ophthalmologic_disease','unenriched','CD31_enriched']:
      z=donor[donor.cells_n.ge(20)].copy()
      if restriction=='no_noted_ophthalmologic_disease':z=z[z.health.eq(restriction)]
      elif restriction in ['unenriched','CD31_enriched']:z=z[z.experiment.eq(restriction)]
      for exclude in z.donor:
        a=z.loc[z.donor.ne(exclude),'positivity_percent']
        lod.append({'restriction':restriction,'omitted_donor':int(exclude),'remaining_donors_n':len(a),
                    'mean_percent':a.mean(),'sample_sd_percent':a.std(ddof=1)})
    pd.DataFrame(lod).to_csv(HERE/'PDE5A_mural_leave_one_donor_out.csv',index=False)
    paired=[]
    for threshold in [20,1,5,10,30,50]:
      z=region[region.cells_n.ge(threshold)].pivot(index='donor',columns='region',values='positivity_percent').dropna()
      for d,r in z.iterrows():
        paired.append({'minimum_cells_each_region':threshold,'donor':int(d),'macula_percent':r['macula'],
                       'peripheral_percent':r['peripheral'],'macula_minus_peripheral_pp':r['macula']-r['peripheral']})
    pd.DataFrame(paired).to_csv(HERE/'PDE5A_mural_paired_region_differences.csv',index=False)
    pd.DataFrame(counts).to_csv(HERE/'sample_cell_count_audit.csv',index=False)
    pd.DataFrame(marker_checks).to_csv(HERE/'PDE5A_and_markers_all_sample_clusters.csv',index=False)
    print(ds[ds.minimum_cells.eq(20)].to_string(index=False))
    return {'total_cells':int(pd.DataFrame(counts).all_cells_n.sum()),'mural_cells':int(donor.cells_n.sum()),
            'annotation':annotation,'primary_summaries':json.loads(ds[ds.minimum_cells.eq(20)].to_json(orient='records'))}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--stage',choices=['markers','quantify'],required=True);args=p.parse_args()
    start=time.perf_counter(); details={}
    if args.stage=='markers':markers()
    else:details=quantify()
    audit={'stage':args.stage,'elapsed_seconds':time.perf_counter()-start,'memory':memory_mb(),
        'python':platform.python_version(),'pandas':pd.__version__,'numpy':np.__version__,
        'protocol_sha256':sha(HERE/'fixed_protocol.json'),'script_sha256':sha(Path(__file__)),
        'source_sha256':{x.name:sha(x) for x in FILES},'details':details}
    (HERE/f'audit_{args.stage}.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
    print('Audit saved',args.stage,round(audit['elapsed_seconds'],2),'seconds',audit['memory'])
