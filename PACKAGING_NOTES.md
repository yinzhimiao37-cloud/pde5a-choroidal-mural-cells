# Archive preparation on 18 November 2026

- Original analysis directory relationships were preserved.
- All analysis rules and archived numerical CSV outputs were retained as copied.
- JSON provenance paths pointing to the author's Windows workspace were made relative. Original source hashes remain in MANIFEST.csv; these text edits do not change archived statistics.
- Two plotting imports were made optional to remove dependence on a personal skill installation. Missing alignment helpers produce an explicit warning. Statistical code was unchanged.
- Two new input-only helpers document exact-expression extraction and hash-checked public Xenium downloads. They were not used to rerun the archived analyses during repository preparation.
- Raw matrices, per-cell spatial/annotation exports, local expression caches, installed runtimes, article PDFs, FAERS record-level inputs, and the historical failed pre-repair LIANA audit were not included. Final LIANA QA and full donor-level contrast tables are retained.
- The new molecular archive is separate from the submitted three-drug FAERS repository; shared reporting data are not described as independent replication.
- Package preparation consisted of file inventory, provenance mapping and static dependency/path review. No new biological analysis or full fresh-environment execution is implied.
- Historical protocol provenance may name an original local metadata-audit path that is not redistributed. The corresponding GEO accession and derived sample metadata are provided; historical protocol text is retained rather than rewritten after analysis.
