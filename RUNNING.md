# Running the molecular modules

Python 3.12 is recommended. Run commands from the repository root after reading
INPUTS.md. Create an isolated environment, then install `requirements.txt`.
For LIANA, use a separate environment if necessary: the archived communication
analysis used pandas 2.3.3, while other archived modules used pandas 3.0.x.
`results/software_versions.json` and module audits preserve original versions.

## Donor localization and cross-dataset consistency

```
python scripts/prepare_pde5a_focus.py
python review_20261013/mural_PDE5A_validation/analyze_GSE183320.py
python review_20261013/mural_PDE5A_validation/GSE135922/run_donor_region_validation.py --stage markers
python review_20261013/mural_PDE5A_validation/GSE135922/run_donor_region_validation.py --stage quantify
python review_20261013/mural_PDE5A_validation/GSE135922/run_metadata_sensitivity.py
```

The saved annotation and protocols are archival inputs. The metadata sensitivity
uses the included `metadata_audit/documented_sample_reuse.csv`.

## Xenium spatial context

```
python scripts/download_spatial_inputs.py
python review_20261013/mural_PDE5A_validation/spatial_GSE318037/preprocess_spatial.py
python review_20261013/mural_PDE5A_validation/analyze_spatial_GSE318037.py
```

The last command runs 9,999 conditional permutations and can take substantial
time. It is not needed to inspect the frozen donor/section/cohort summaries.

## State panels and diffusion pseudotime

```
python review_20261013/mural_heterogeneity_extension/scRNA/prepare_expression.py
python review_20261013/mural_heterogeneity_extension/scRNA/analyze_mural_states.py
```

The generated AnnData cache is also the default input for LIANA. It is not
redistributed. The expression preparation preserves published normalization.

## LIANA donor-paired expression support

```
python review_20261013/mural_heterogeneity_expansion/communication/run_liana_donor.py --stage all
python review_20261013/mural_heterogeneity_expansion/communication/audit_liana_values.py
```

Review the archived protocol for group definitions, donor thresholds, matching,
multiple-testing families and unevaluable contrasts. Do not select only
unmatched positive contrasts from the full results.

## Plotting from included summary tables

```
python review_20261013/mural_PDE5A_validation/plot_validation.py
python review_20261013/mural_PDE5A_validation/plot_spatial.py
```

These commands overwrite figure exports. The author's optional personal
panel-alignment helper is no longer mandatory. If absent, a warning is emitted
and that display audit is skipped; no claim of a new alignment audit is made.
The archived figure PDFs remain available for direct inspection.
