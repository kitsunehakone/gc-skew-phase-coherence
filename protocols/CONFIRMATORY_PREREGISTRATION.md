# Confirmatory Cohort A — Frozen preregistration

## Purpose
Independent confirmation of the higher-odd phase-coherence signal discovered in the opened nine-genome pilot.

The opened nine-genome pilot and the older 15-genome development set are not reused as confirmation data.

## Cohort
Primary cohort: 15 chromosomes, one per genus.

Eligibility was frozen before downloading or inspecting sequence content:
1. bacterial chromosome listed by DoriC;
2. complete assembly;
3. circular topology;
4. a single DoriC-listed replication-origin region on the chromosome;
5. genus absent from both prior opened sets;
6. one genome per genus;
7. exact chromosome accession recorded before analysis.

Five reserve chromosomes are frozen separately. A reserve may replace a primary genome only for a provenance/eligibility failure, never because of a result.

## Sequence gate
No GC-skew statistic is computed until:
- FASTA parses successfully;
- FASTA header matches the frozen accession;
- sequence contains only IUPAC DNA symbols;
- sequence length equals the frozen DoriC chromosome length;
- raw and canonical SHA-256 hashes are recorded.

## Frozen analysis
Primary resolution: 512 bins.
Robustness resolutions: 256 and 1024 bins.

Per chromosome:
- J_max
- first-harmonic energy P1
- half-turn-odd energy E_Hminus
- higher-odd energy E_Hminus - P1
- odd-axis coherence C_odd = J_max/E_Hminus
- four (R,H) sector fractions
- symmetry-sector entropy

### Primary null
B = 4,999 Fourier surrogates per chromosome.
Every surrogate preserves:
- every Fourier magnitude exactly;
- the full complex k=1 coefficient exactly;
and randomizes only higher-harmonic phases.

p_i = (1 + #{J_null >= J_obs})/(B+1).

### Primary cohort-level endpoint
Fisher combination:
X = -2 sum_i log(p_i), compared with chi-square on 30 degrees of freedom.

The confirmatory study is positive iff:
p_Fisher < 0.05.

This rule is frozen before any Confirmatory Cohort A sequence is opened.

## Secondary endpoints
- number/proportion with raw p_i < 0.05;
- Wilson 95% CI for that proportion;
- BH q-values across 15 tests;
- number with q < 0.05;
- higher-odd energy and C_odd distributions;
- 256/512/1024 resolution stability.

Localization superiority is not a primary endpoint.

## No post-opening tuning
After any primary sequence is opened, do not alter:
- J;
- B;
- primary resolution;
- the k=1-preserving null;
- cohort membership except by the frozen reserve rule;
- the Fisher p<0.05 success criterion.
