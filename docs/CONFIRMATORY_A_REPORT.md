# Confirmatory Cohort A — Frozen phase-coherence result

## Preregistered outcome

The confirmatory cohort contained **15** previously unopened complete circular
bacterial chromosomes, one per genus. Every chromosome passed the frozen
provenance gate before GC-skew analysis.

The primary null used **4,999** Fourier surrogates per chromosome. Each
surrogate preserved every Fourier magnitude and the complete complex
first-harmonic coefficient while randomizing only higher-harmonic phases.

The prespecified cohort-level test was Fisher's combination of the 15
one-sided Monte Carlo p-values.

## Primary result

$$
X = 104.456328, \qquad df = 30,
$$

with

$$
\boxed{p_{\rm Fisher} = 3.62659e-10}.
$$

The frozen success criterion was $p_{\rm Fisher}<0.05$.

Therefore:

$$
\boxed{\text{Confirmatory Cohort A is POSITIVE for higher-harmonic phase coherence.}}
$$

This conclusion uses the criterion specified before any confirmatory sequence
was opened.

## Prevalence

**6/15 = 40.0%**
chromosomes had raw genome-level $p<0.05$.

The Wilson 95% confidence interval for this proportion is

$$
[0.198,\,0.643].
$$

After Benjamini-Hochberg correction, **6/15**
chromosomes remained significant at $q<0.05$.

## BH-significant chromosomes

- **Streptococcus mitis NCTC 12261** (`NZ_CP028414.1`): $p=0.0002$, $q=0.003$, $J_{\max}=0.8425$, $C_{\rm odd}=0.9569$, higher-odd energy=0.2604, resolution-stable=yes.
- **Enterococcus casseliflavus EC20** (`NC_020995.1`): $p=0.0004$, $q=0.003$, $J_{\max}=0.8738$, $C_{\rm odd}=0.9554$, higher-odd energy=0.2191, resolution-stable=yes.
- **Acinetobacter baumannii D1279779** (`NC_020547.2`): $p=0.0006$, $q=0.003$, $J_{\max}=0.7863$, $C_{\rm odd}=0.9308$, higher-odd energy=0.2515, resolution-stable=yes.
- **Clostridioides difficile CF5** (`NC_017173.1`): $p=0.0008$, $q=0.003$, $J_{\max}=0.8655$, $C_{\rm odd}=0.9613$, higher-odd energy=0.2045, resolution-stable=no.
- **Thioalkalivibrio paradoxus ARh 1** (`NZ_CP007029.1`): $p=0.0038$, $q=0.0114$, $J_{\max}=0.7031$, $C_{\rm odd}=0.8980$, higher-odd energy=0.2827, resolution-stable=yes.
- **Cytobacillus oceanisediminis 2691** (`NZ_CP015506.1`): $p=0.0054$, $q=0.0135$, $J_{\max}=0.8519$, $C_{\rm odd}=0.9609$, higher-odd energy=0.1741, resolution-stable=yes.

## Cohort summaries

- Median $J_{\max}$: **0.6236**
- Median first-harmonic energy $P_1$: **0.4968**
- Median higher-odd energy fraction: **0.2827**
- Median odd-axis coherence $C_{\rm odd}$: **0.8859**
- Resolution-stable under the frozen 256/512/1024 criterion: **10/15**

## Interpretation

The result supports the narrow preregistered claim that, across this independent
cohort, bacterial GC-skew profiles contain higher-harmonic phase organization
that is not explained by the Fourier magnitude spectrum or by the complete
$k=1$ coefficient alone.

This is stronger than observing higher odd-frequency power. The null preserves
that power and destroys only higher-harmonic phase relations.

The result does **not** establish that the projector universally improves oriC
localization, does not identify a new molecular mechanism, and does not imply a
thermodynamic or symbolic interpretation.

The cohort contains one chromosome per genus to reduce obvious duplication,
but genomes are not statistically independent evolutionary replicates in the
strict phylogenetic sense. The Fisher result should therefore be accompanied
by the per-genome p-values, prevalence estimate, and a later phylogenetically
broader analysis.

## Tomography endpoint

The separately frozen tomography statistic $s_4/s_1$ has **not** been computed
here. The exact pre-existing tomography operator/fiber construction is not
present in the currently available executable files. Because the confirmatory
sequences are now opened, that construction must not be invented or altered
post hoc. It should be run only from the already-defined tomography
implementation/specification that existed before sequence opening.
