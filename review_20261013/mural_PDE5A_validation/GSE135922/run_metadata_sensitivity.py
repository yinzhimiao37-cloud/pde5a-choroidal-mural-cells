"""Metadata-driven additional sensitivity; does not alter the fixed primary analysis.

Added after the primary run when sample reuse and cross-source clinical-label
conflicts were documented. These exclusions are not selected using PDE5A.
"""
import hashlib,json
from pathlib import Path
import pandas as pd
HERE=Path(__file__).resolve().parent
source=HERE/'PDE5A_mural_donor_pooled_regions.csv'
reuse_file=HERE.parent/'metadata_audit/documented_sample_reuse.csv'
d=pd.read_csv(source)
r=pd.read_csv(reuse_file)
shared=set(r.loc[r.source_dataset.eq('GSE135922'),'GSE135922_donor_label'].astype(int))
assert shared=={4,5,6,7}, shared
rules={
 'exclude_documented_shared_donors':~d.donor.isin(shared),
 'original_no_noted_disease_exclude_discordant_donors5_6':d.health.eq('no_noted_ophthalmologic_disease')&~d.donor.isin([5,6])
}
points=[];stats=[];loo=[]
for label,mask in rules.items():
 x=d[mask].copy();x['restriction']=label;points.append(x)
 for threshold in [20,1,5,10,30,50]:
  z=x[x.cells_n.ge(threshold)];a=z.positivity_percent
  stats.append({'restriction':label,'minimum_cells':threshold,'donors_n':len(z),'donor_ids':';'.join(map(str,z.donor)),
    'cells_n':int(z.cells_n.sum()),'positive_n':int(z.positive_n.sum()),'mean_percent':a.mean(),
    'sample_sd_percent':a.std(ddof=1),'donors_with_positive_cells':int(z.positive_n.gt(0).sum())})
 for omit in x[x.cells_n.ge(20)].donor:
  a=x.loc[x.cells_n.ge(20)&x.donor.ne(omit),'positivity_percent']
  loo.append({'restriction':label,'omitted_donor':int(omit),'remaining_donors_n':len(a),
               'mean_percent':a.mean(),'sample_sd_percent':a.std(ddof=1)})
pd.concat(points).to_csv(HERE/'metadata_sensitivity_donor_points.csv',index=False)
pd.DataFrame(stats).to_csv(HERE/'metadata_sensitivity_summary.csv',index=False)
pd.DataFrame(loo).to_csv(HERE/'metadata_sensitivity_leave_one_out.csv',index=False)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
audit={'timing':'added after primary statistics on discovery of accession-level sample reuse and clinical-label discordance',
       'selection_basis':'metadata only, not PDE5A results',
       'documented_shared_original_donors':sorted(shared),
       'source_specific_no_noted_disease_map':'retained in primary analysis; not asserted healthy across all publications',
       'conservative_health_exclusion':'exclude both original donor5 and donor6 because cross-paper age/sex and revised AMD labels cannot be reconciled confidently',
       'input_sha256':{p.name:sha(p) for p in [source,reuse_file]},'script_sha256':sha(Path(__file__))}
(HERE/'audit_metadata_sensitivity.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
print(pd.DataFrame(stats).query('minimum_cells==20').to_string(index=False))
