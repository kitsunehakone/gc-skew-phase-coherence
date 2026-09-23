# Provenance timeline

This repository intentionally contains materials from both sides of the confirmatory unblinding boundary.

## Fixed before Confirmatory Cohort A sequence opening

- `code/gc_symmetry_projector.py`
- `code/skewdb_projector_benchmark.py`
- `code/run_confirmatory_prevalence.py`
- `protocols/CONFIRMATORY_PREREGISTRATION.md`
- `protocols/TOMOGRAPHY_STATISTIC_ADDENDUM_v1_0.md`
- the primary and reserve accession manifests contained in the frozen confirmatory bundle

The three files in `code/` are byte-for-byte identical to the copies in the pre-unblinding frozen bundle. Their SHA-256 values are listed in `FROZEN_CODE_SHA256SUMS`.

## Produced after sequence opening

- `results/`
- `figures/`
- confirmatory provenance output under `manifests/confirmatory_A_provenance_gate.csv`
- manuscript/report material in `docs/`
- `tools/run_confirmatory_chunked.py` (execution convenience wrapper only)

The post-unblinding chunked wrapper does not define or modify the statistic, null, cohort, number of surrogates, seed schedule, or success criterion. It invokes the frozen runner one accession at a time and recomputes only the predeclared cohort summaries after all jobs complete.

## Unexecuted prespecified endpoint

The tomography `s4/s1` endpoint was frozen before sequence opening but was not executed because its pre-existing implementation was unavailable at unblinding. It was not reconstructed after the cohort was opened.
