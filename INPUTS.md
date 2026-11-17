# Input data and provenance

## GSE183320

GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183320

Download the author metadata and full normalized expression matrix into:

```
data/raw/GSE183320_all_metadata.csv.gz
data/raw/GSE183320_all_normalized_expression_data.csv.gz
```

The localization script also uses `data/processed/sc_choroid_focus_matrix.csv`,
with columns `cell_barcode` and `PDE5A`. Generate this exact-expression extract
using `python scripts/prepare_pde5a_focus.py` after downloading the two inputs.
This helper extracts the published PDE5A row without normalization or imputation.
The original expression and metadata hashes are recorded in
`review_20261013/mural_PDE5A_validation/GSE183320/analysis_audit.json`.
The original focused-matrix hash may differ if the helper exports a different
floating-point text representation; the analysis directly checks the extracted
values against the original full matrix before analysis.

## GSE135922

GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135922

Place all 14 `GSM4037*expression.tsv.gz` files in `data/raw/` with their original
names. Filenames and SHA256 values are recorded in `GSE135922/audit_quantify.json`.
`fixed_protocol.json`, `annotation.json`, and the metadata-overlap evidence table
are included. The marker stage summarizes lineages; it does not authorize
relabeling based on PDE5A positivity.

## GSE318037 Xenium

GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318037

The archived `spatial_GSE318037/input/download_manifest.json` lists 22 public
URLs, filenames, sizes and hashes: seven feature matrices, seven cell-summary
files, seven feature lists, and the author merged annotation. Run
`python scripts/download_spatial_inputs.py` to retrieve exactly these files into
the expected input directory. This helper validates both bytes and SHA256 and
does not modify the archived manifest. Then run `preprocess_spatial.py` to
create the joined cell table locally. Coordinates are centroids in micrometres.
Seven sections represent six donors; CNV3A and CNV3B are one donor.

## FAERS

Official raw quarterly releases must be obtained from the FDA. The focused
aggregate scenario tables included here are not substitutes for a report-level
linkage audit. Shared three-drug extraction code and its prerequisites are in
the separate FAERS repository linked in README.md.
