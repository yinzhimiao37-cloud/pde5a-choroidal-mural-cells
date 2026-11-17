"""Join public Xenium coordinates, PDE5A counts and published cell annotation.

No expression or neighborhood-based inclusion is performed here. All raw
segmented cells remain in the output; unavailable published annotation is flagged.
"""
from pathlib import Path
import gzip, io, json, hashlib
import h5py, numpy as np, pandas as pd
from scipy.sparse import csc_matrix

OUT=Path(__file__).resolve().parent
INPUT=OUT/'input'
SAMPLES=[
 ('GSM9485163','control_donor_1','Control_donor_1','Control1','Control1','non-CNV'),
 ('GSM9485164','control_donor_2','Control_donor_2','Control2','Control2','non-CNV'),
 ('GSM9485165','control_donor_3','Control_donor_3','Control3','Control3','non-CNV'),
 ('GSM9485166','CNV_donor_1','CNV_donor_1','CNV1','CNV1','CNV'),
 ('GSM9485167','CNV_donor_2','CNV_donor_2','CNV2','CNV2','CNV'),
 ('GSM9485168','CNV_donor_3A','CNV_donor_3A','CNV3A','CNV3','CNV'),
 ('GSM9485169','CNV_donor_3B','CNV_donor_3B','CNV3B','CNV3','CNV'),
]

def decode(arr):return np.array([x.decode() if isinstance(x,bytes) else str(x) for x in arr])

def main():
    meta=pd.read_csv(INPUT/'GSE318037_merged_xenium_metadata.csv.gz')
    if meta.barcodes.duplicated().any():raise ValueError('Published metadata barcode duplicates')
    result=[];audit=[];panels=[]
    for gsm,prefix,original,section,donor,condition in SAMPLES:
        with gzip.open(INPUT/f'{gsm}_{prefix}_cells.parquet.gz','rb') as f:
            cells=pd.read_parquet(io.BytesIO(f.read()))
        if cells.cell_id.duplicated().any():raise ValueError('Duplicated raw cells')
        with h5py.File(INPUT/f'{gsm}_{prefix}_cell_feature_matrix.h5','r') as h:
            group=h['matrix'];names=decode(group['features/name'][:]);ids=decode(group['features/id'][:])
            types=decode(group['features/feature_type'][:]);barcodes=decode(group['barcodes'][:])
            ix=np.flatnonzero((names=='PDE5A') & (ids=='ENSG00000138735') & (types=='Gene Expression'))
            if len(ix)!=1:raise ValueError(f'{gsm}: PDE5A feature not uniquely defined')
            matrix=csc_matrix((group['data'][:],group['indices'][:],group['indptr'][:]),shape=tuple(group['shape'][:]))
            vals=np.asarray(matrix[ix[0],:].toarray()).ravel()
            gene_counts=np.asarray(matrix[types=='Gene Expression'].sum(axis=0)).ravel()
        counts=pd.DataFrame({'cell_id':barcodes,'PDE5A_count':vals,'matrix_gene_transcripts':gene_counts})
        if counts.cell_id.duplicated().any():raise ValueError('Duplicated matrix cell ids')
        if set(counts.cell_id)!=set(cells.cell_id):raise ValueError('Matrix and cells ID sets differ')
        joined=cells.merge(counts,on='cell_id',how='left',validate='one_to_one')
        joined['barcode']=original+'_'+joined.cell_id
        sample_meta=meta[meta.donor_id.eq(original)].copy()
        sample_meta=sample_meta.rename(columns={'donor_id':'source_section','barcodes':'barcode',
                                                'final_celltype':'celltype','tissue':'region',
                                                'segmentation_method':'annotation_segmentation_method'})
        if not set(sample_meta.barcode).issubset(joined.barcode):raise ValueError('Annotation ID not found in coordinate table')
        joined=joined.merge(sample_meta,on='barcode',how='left',validate='one_to_one',indicator=True)
        joined['author_annotation_available']=joined['_merge'].eq('both');joined=joined.drop(columns='_merge')
        joined['section_id']=section;joined['donor_id']=donor;joined['condition']=condition;joined['gsm']=gsm
        joined=joined.rename(columns={'x_centroid':'x_um','y_centroid':'y_um',
                                      'transcript_counts':'total_transcripts','total_counts':'total_codewords'})
        # Area/QC flag is a transparent data-validity check, not a candidate-selection filter.
        joined['coordinate_qc']=np.isfinite(joined.x_um)&np.isfinite(joined.y_um)
        joined['area_qc']=np.isfinite(joined.cell_area)&joined.cell_area.gt(0)&np.isfinite(joined.nucleus_area)&joined.nucleus_area.ge(0)
        joined['nucleus_within_cell_area']=joined.nucleus_area.le(joined.cell_area)
        joined['single_nucleus']=joined.nucleus_count.eq(1)
        joined['celltype']=joined.celltype.fillna('unannotated');joined['region']=joined.region.fillna('unannotated')
        result.append(joined)
        audit.append({'section_id':section,'donor_id':donor,'gsm':gsm,'raw_cells':len(cells),
                      'matrix_cells':len(counts),'annotated_cells':len(sample_meta),
                      'unannotated_raw_cells':int((~joined.author_annotation_available).sum()),
                      'all_metadata_ids_found':True,
                      'matrix_gene_sum_vs_cell_transcripts_exact_matches':int(joined.matrix_gene_transcripts.eq(joined.total_transcripts).sum()),
                      'coordinate_qc_pass':int(joined.coordinate_qc.sum()),'area_qc_pass':int(joined.area_qc.sum())})
        panels.append({'gsm':gsm,'total_features':len(names),'gene_expression_features':int((types=='Gene Expression').sum()),
                       'PDE5A_feature_index_zero_based':int(ix[0]),'PDE5A_feature_id':ids[ix[0]],
                       'feature_names_sha256':hashlib.sha256('\n'.join(names).encode()).hexdigest()})
    allcells=pd.concat(result,ignore_index=True)
    if allcells.duplicated(['section_id','cell_id']).any():raise ValueError('Nonunique section and cell IDs')
    primary=['section_id','donor_id','condition','gsm','cell_id','barcode','region','celltype','x_um','y_um',
             'PDE5A_count','total_transcripts','cell_area','nucleus_area','nucleus_count','area_qc',
             'coordinate_qc','author_annotation_available','nucleus_within_cell_area','single_nucleus',
             'segmentation_method','total_codewords','matrix_gene_transcripts','nCount_Xenium','nFeature_Xenium']
    allcells=allcells[primary+[x for x in allcells if x not in primary]]
    allcells.to_parquet(OUT/'joined_all_cells.parquet',index=False)
    allcells.to_csv(OUT/'joined_all_cells.csv.gz',index=False)
    pd.DataFrame(audit).to_csv(OUT/'join_audit_by_section.csv',index=False)
    pd.DataFrame(panels).to_csv(OUT/'PDE5A_feature_coverage.csv',index=False)
    allcells.groupby(['section_id','donor_id','region','celltype'],observed=True).size().rename('n_cells').reset_index().to_csv(OUT/'region_celltype_counts.csv',index=False)
    summary={'raw_cells':len(allcells),'annotated_cells':int(allcells.author_annotation_available.sum()),
             'sections':7,'biological_donors':6,'all_output_columns':list(allcells),
             'celltype_labels':sorted(allcells.celltype.unique()),'region_labels':sorted(allcells.region.unique()),
             'total_transcripts_definition':'cells.parquet transcript_counts: gene-expression transcripts; total_codewords also includes controls',
             'PDE5A_definition':'Unique ENSG00000138735 / PDE5A Gene Expression row in raw per-cell feature matrix; no normalization',
             'annotation_origin':'Original authors label transfer from earlier scRNA-seq; tissue region manually traced in Xenium Explorer',
             'area_qc_definition':'finite cell_area>0 and finite nucleus_area>=0; no row removed',
             'no_independent_cell_annotation_validation':True,
             'coordinate_units':'Micrometers under Xenium cells.parquet standard; evidence documented separately',
             'no_spatial_association_test_in_preprocessing':True}
    (OUT/'preprocessing_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps({x:summary[x] for x in ['raw_cells','annotated_cells','sections','biological_donors','celltype_labels','region_labels']},indent=2))

if __name__=='__main__':main()
