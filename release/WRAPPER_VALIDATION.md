# Chunked-wrapper validation

The post-unblinding convenience wrapper was validated against the archived confirmatory output at the full `4,999`-surrogate setting on two chromosomes spanning the cohort order and genome-size range:

- `NC_017463.1` (cohort row 1; adjusted seed `20260920`)
- `NZ_CP015506.1` (cohort row 15; adjusted seed `21660962`)

For both accessions, the wrapper reproduced the archived values exactly (floating-point difference `0.0`) for:

- `J_max`
- `k1_preserved_phase_p`
- `projector_axis_fraction`
- `first_harmonic_power_fraction`
- `higher_odd_energy_fraction`
- `odd_axis_coherence`

The wrapper's cohort-statistic helper functions were also checked against the archived 15-row table:

- Fisher statistic: `104.45632769223512`
- Fisher survival probability: `3.6265898370043044e-10`
- Benjamini-Hochberg q-values: exact numerical match to archived values.

This validation concerns execution equivalence only. The wrapper remains post-unblinding infrastructure and is not part of the prespecified scientific analysis.

## v3 merge-integrity guards

The post-unblinding wrapper was additionally checked after repository v2:

- `chi2_sf_even` now rejects odd or non-positive degrees of freedom; the cohort path still uses the exact even-df closed form (`df=30`).
- A valid one-genome merge completed successfully from an archived full-precision row.
- A deliberately truncated result row with the required `odd_axis_coherence` field removed failed loudly with `Truncated/incomplete result`, as intended.
- Merge also now requires a one-row provenance file, exact accession agreement, finite required numerical fields, and a Monte Carlo p-value within `[1/(B+1), 1]`.

These changes are wrapper-only integrity checks and do not modify the frozen pre-unblinding scientific runner.
