# Repository release v3

Repository-only reproducibility hardening. No frozen scientific code, statistic, null model, cohort, endpoint, result, or manuscript claim was changed.

Changes from repository package v2:

- clarified that the full 15-genome 4,999-surrogate run is expected to finish in roughly under three minutes on comparable hardware; the prior interrupted monolithic rerun reflected an interactive wall-clock limit, not computational heaviness;
- hardened the post-unblinding chunked merge wrapper to fail loudly on missing provenance, missing required result columns, accession mismatches, non-finite required values, or impossible Monte Carlo p-values;
- added an explicit positive-even-df guard to the wrapper-only `chi2_sf_even` helper;
- removed generated Python `__pycache__` files from the public package;
- regenerated release-wide checksums.

The frozen pre-unblinding files in `code/` and `protocols/` remain byte-for-byte unchanged from v2 and from the archived pre-unblinding bundle.
