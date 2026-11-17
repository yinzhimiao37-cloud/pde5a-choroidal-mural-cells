"""结果后仅作配体家族注释；不重新筛选/校正，也不把LR的q叫通路q。"""
from pathlib import Path
import pandas as pd
HERE=Path(__file__).resolve().parent
s=pd.read_csv(HERE/'results/donor_LR_statistics_all_families.csv')
families={
 'TGF-beta ligands':{'TGFB1','TGFB2','TGFB3'},
 'PDGF ligands':{'PDGFA','PDGFB','PDGFC','PDGFD'},
 'Notch ligands':{'JAG1','JAG2','DLL1','DLL3','DLL4'},
 'Endothelin ligands':{'EDN1','EDN2','EDN3'},
 'Angiopoietin ligands':{'ANGPT1','ANGPT2','ANGPT4'},
 'VEGF ligands':{'VEGFA','VEGFB','VEGFC','FIGF'}}
rows=[]
for name,genes in families.items():
 t=s[s.ligand_complex.map(lambda x: bool(set(str(x).split('_'))&genes))].copy()
 t['ligand_family_annotation']=name
 rows.append(t)
out=pd.concat(rows,ignore_index=True)
out.to_csv(HERE/'results/preidentified_ligand_family_lookup.csv',index=False)
summary=out.groupby(['ligand_family_annotation','variant','cohort']).agg(
   evaluable_LR=('q_BH','size'),positive_candidates=('candidate_support','sum'),
   smallest_LR_q_not_a_pathway_q=('q_BH','min')).reset_index()
summary.to_csv(HERE/'results/ligand_family_descriptive_summary.csv',index=False)
print('Annotated',len(out),'rows; original full-family q values retained.')
