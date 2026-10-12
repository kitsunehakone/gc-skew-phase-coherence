# Blinded Primitive-Closure Reserve Validation — Frozen Protocol

Frozen before downloading or inspecting any reserve-sequence FASTA.

## Purpose
Independent blinded validation of the newly derived primitive harmonic-closure prediction on the five reserve chromosomes that were frozen before Confirmatory Cohort A was opened and were not used in the confirmatory analysis.

## Cohort
Exactly the five accessions listed in `manifests/CONFIRMATORY_RESERVE_QUEUE.csv`.
No replacement, exclusion, reordering, or substitution is permitted for statistical reasons.

## Sequence gate
Before any GC-skew statistic is accepted:
1. FASTA must parse successfully.
2. Header/accession must correspond to the frozen accession.
3. Sequence must contain only IUPAC DNA symbols.
4. Sequence length must equal the frozen expected chromosome length.
5. Raw and canonical sequence SHA-256 hashes are recorded.

A provenance failure invalidates that genome; it is not replaced.

## Signal
Primary resolution M=512 equal contiguous bins.
Within each bin:
  s=(G-C)/(G+C), with zero only if G+C=0.
Mean-center s before Fourier analysis.

Robustness resolutions M=256 and M=1024 are secondary only.

## Frozen primitive-closure statistic
Let c_k be the positive-frequency Fourier coefficients of the centered real GC-skew signal.
For odd n=3,5,... with n+2 below Nyquist, define

  kappa_n = arg(c_{n+2}) + arg(c_{n-2}) - 2 arg(c_n)

and weights

  w_n = |c_{n-2}| |c_n|^2 |c_{n+2}|.

Define

  M1 = [sum_n w_n exp(i kappa_n)] / [sum_n w_n]

and the directional primitive-closure statistic

  S1 = Re(M1).

The directional alternative is S1 larger than expected under the null, corresponding to primitive closure phases preferentially aligned near kappa=0.

No statistic, weight, harmonic range, direction, or threshold may be altered after reserve FASTAs are opened.

## Frozen null
B=4,999 surrogates per genome.

Every surrogate preserves:
- every Fourier magnitude exactly;
- the full complex k=1 coefficient exactly;
- real-signal conjugate symmetry;

and randomizes all higher positive-frequency phases.

Per-genome one-sided Monte Carlo p-value:

  p_i = [1 + #{S1_null >= S1_obs}] / 5000.

Minimum attainable p is 0.0002.

## Primary cohort-level endpoint
Combine the five per-genome p-values by Fisher's method:

  X = -2 sum_i log p_i

and compare with chi-square on 10 degrees of freedom.

The blinded reserve validation is declared positive iff

  p_Fisher < 0.05.

This is the sole primary success criterion.

## Secondary / descriptive outputs
- per-genome S1 and p_i;
- Benjamini-Hochberg q-values across the five tests;
- number with raw p<0.05 and q<0.05;
- |M1|, arg(M1), effective closure-weight count;
- sign-quotiented M2 for descriptive comparison only;
- M=256 and M=1024 repetitions as robustness checks;
- reproduction/provenance diagnostics.

Neither R1=|M1| nor M2 may replace S1 as the primary endpoint after unblinding.

## Interpretation
This is an independent test of the primitive closure prediction on previously unused reserve chromosomes. It does not establish microscopic equivalence between bacterial replication and harmonic light. A positive result would independently validate the predicted GC-side primitive closure alignment under the same k=1-preserving, full-magnitude-preserving null.
