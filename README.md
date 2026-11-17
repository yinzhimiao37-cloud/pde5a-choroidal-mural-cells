# PDE5A heterogeneity in human choroidal mural cells

Code, frozen analysis rules, donor-level summary data, and publication figures for
the manuscript **PDE5A heterogeneity in human choroidal mural cells and
source-dependent tadalafil reporting in FAERS**, prepared for Graefe's Archive
for Clinical and Experimental Ophthalmology. This description does not imply acceptance.

Repository: https://github.com/yinzhimiao37-cloud/pde5a-choroidal-mural-cells

## What is included

| Module | Directory | Principal entry point |
|---|---|---|
| GSE183320 donor localization and sampling/depth sensitivity | `review_20261013/mural_PDE5A_validation/` | `analyze_GSE183320.py` |
| GSE135922 donor and regional consistency | `review_20261013/mural_PDE5A_validation/GSE135922/` | `run_donor_region_validation.py` |
| GSE318037 Xenium preprocessing and conditional spatial permutations | `review_20261013/mural_PDE5A_validation/` | `analyze_spatial_GSE318037.py` |
| Mural-cell clustering, expression panels and diffusion pseudotime | `review_20261013/mural_heterogeneity_extension/scRNA/` | `prepare_expression.py`, `analyze_mural_states.py` |
| Donor-paired LIANA expression-support analysis | `review_20261013/mural_heterogeneity_expansion/communication/` | `run_liana_donor.py` |

The original relative directory layout is retained because the analysis scripts
resolve inputs and outputs relative to their own location. Summary CSV/JSON files
are stored beside their original modules. Four main-figure PDFs are in
`review_20261013/mural_PDE5A_validation/figures/`.

## Reproducibility and input access

This is a code and summary-data archive, not a zero-setup one-click pipeline.
Published normalized expression values are not raw UMI counts. Full expression
matrices, per-cell coordinate tables, expression caches, raw FAERS releases,
report identifiers and local databases are not redistributed here.
See [INPUTS.md](INPUTS.md) for source accessions, required paths and preparation,
and [RUNNING.md](RUNNING.md) for module-specific commands. Archived audits record
the versions and hashes used for the original results. The general requirements
file is an installation starting point, not an assertion that all future package
versions will reproduce identical graphs. For archival versions consult each
module's JSON audit and original requirements file.

Viewing the included summary tables needs no raw data. Running an analysis
writes outputs into that module; use a separate clone or copy to preserve the
frozen tables before rerunning. Repository preparation did not rerun the large
analyses or regenerate their numerical results.

## FAERS context and related manuscript

The molecular study and a related three-drug ocular FAERS manuscript use
overlapping FAERS releases and processing. They are not independent FAERS
replications. Shared FAERS spectrum code is archived separately at
https://github.com/yinzhimiao37-cloud/pde5i-ocular-faers.
The four `Audit_Data_*.csv` files under
`review_20261009/revision_draft6/source_data/` contain aggregate focused
tadalafil-CSCR mapping and source-sensitivity scenarios. They permit inspection
of the supplied counts, not reconstruction of patient identity or a complete
independent reproduction of all endpoint-linkage decisions. This repository
does not include the full focused-FAERS case-level evidence ledger, original
record-level linkage inputs, or all full-background extraction code.

## Interpretation

- Biological replication is at donor level; effective donor counts differ by population and eligibility rule.
- GSE135922 and GSE183320 reuse four documented donors. Non-shared subsets are explicitly identified.
- Spatial material is non-CSCR tissue without established tadalafil exposure. A possible spatial Control3 overlap remains unresolved.
- Frozen primary spatial analyses did not support positive endothelial-neighbour enrichment.
- Pseudotime does not establish chronological ageing. Transcript-based ligand–receptor scores are not measured cellular responses; matching removed adjusted support for the positive contrasts described in the manuscript.
- PDE5A transcript detection does not measure protein activity or reproduce tadalafil inhibition.
- Analyses were post hoc extensions; dated protocols document what was fixed before new summaries, not preregistration before the original exploration.

## Packaging notes and licensing

See [PACKAGING_NOTES.md](PACKAGING_NOTES.md) and [MANIFEST.csv](MANIFEST.csv).
Original analysis source hashes and packaged file hashes distinguish copies
from display-only portability edits. Archived numerical results were copied,
not recalculated. Author code is provided under the MIT license. Public source
datasets and third-party interaction resources retain their original terms;
this license does not grant new rights over those sources. Cite the original
dataset publications and LIANA when reusing their resources.
