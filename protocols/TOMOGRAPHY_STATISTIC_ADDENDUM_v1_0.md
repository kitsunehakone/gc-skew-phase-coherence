# Confirmatory Cohort A — Tomography Statistic Addendum v1.0

## Status

This addendum is frozen **before any Confirmatory Cohort A sequence is opened**.
It does not modify the previously frozen cohort membership, reserve queue,
projector statistic J, phase-null construction, or Fisher cohort-level rule.

It adds one prespecified tomography endpoint.

## Frozen tomography construction

For each genome, construct the tomography matrix/operator using the already
defined tomography pipeline with no post-opening changes to:

- binning,
- antipodal pairing,
- fiber construction,
- even/odd channel definitions,
- normalization,
- matrix dimensions,
- or singular-value computation.

Let its singular values, in nonincreasing order, be

    s1 >= s2 >= s3 >= s4 >= ...

## Primary tomography statistic

The primary tomography statistic is

    T = s4 / s1.

Lower values indicate a sharper spectral break after the first three singular
directions.

This ratio is chosen instead of an algorithmically detected "knee location"
because it is continuous, dimensionless, uniquely defined once the tomography
matrix is fixed, and does not introduce a tunable knee-detection rule.

No alternate ratio (for example s3/s1, s5/s1, s4/s2) may replace T after the
confirmatory sequences are opened.

## Primary phase-randomized null

For each genome, generate B = 4,999 surrogate GC-skew signals.

Each surrogate preserves:

1. every Fourier magnitude exactly;
2. the complete complex k=1 coefficient exactly;
3. conjugate symmetry required for a real signal.

Only the phases of higher harmonics are randomized.

For every surrogate, recompute the complete tomography construction and its
singular values, then calculate

    T_null = s4_null / s1_null.

The one-sided per-genome p-value is

    p_i = (1 + #{T_null <= T_obs}) / (B + 1).

The alternative is therefore prespecified as an unusually *small* s4/s1 ratio,
i.e. a sharper rank break than expected from the frozen spectrum and k=1
component alone.

## Cohort-level confirmatory endpoint

For the 15 frozen primary genomes, combine the 15 one-sided p-values using the
already frozen Fisher rule:

    X = -2 * sum_i log(p_i).

The tomography confirmation criterion is

    p_Fisher < 0.05.

This is a separate confirmatory endpoint from the previously frozen J-based
phase-coherence endpoint. Both outcomes must be reported regardless of whether
either is significant.

## Secondary reporting

Report, without changing the primary criterion:

- all 15 observed T = s4/s1 values;
- each genome's null median and 95% null interval for T;
- raw p_i;
- Benjamini-Hochberg q-values across the 15 tomography tests;
- the number and proportion of genomes with raw p_i < 0.05;
- the number with BH q < 0.05;
- resolution robustness at 256, 512, and 1024 bins.

## Knee location

"Knee location" is not a confirmatory endpoint in Cohort A.

It may be shown descriptively in figures after the primary analysis, but no
knee-based hypothesis test or success claim may replace or override the frozen
T = s4/s1 result.

## Interpretation

A significant result would support the specific statement that the observed
tomography spectra exhibit a sharper fourth-singular-value collapse than is
expected from the Fourier power spectrum and the complete first harmonic alone.

It would not, by itself, establish a new biological mechanism, thermodynamic
claim, symbolic correspondence, or universal law.

## No symbolic mapping

The previously discussed hand-drawn symbols are excluded from the mathematical
and empirical argument. No symbol-to-sector correspondence is part of the
theorem, statistic, null model, or confirmatory analysis.
