# Symmetry-resolved GC-skew phase coherence

Reproducibility package for the manuscript **“Symmetry-resolved GC-skew reveals higher-harmonic phase coherence beyond the first Fourier mode in bacterial chromosomes.”**

## What this repository supports

The confirmatory analysis tests whether higher odd Fourier harmonics of bacterial GC-skew profiles show common-axis phase organization after conditioning on:

1. the complete Fourier magnitude spectrum, and
2. the exact complex first-harmonic (`k=1`) coefficient.

The primary statistic is the maximum reversal-odd / half-turn-odd projector energy fraction, `J_max`. The predeclared confirmatory cohort consists of 15 complete circular bacterial chromosomes from 15 genera. The frozen cohort-level endpoint was Fisher's combination of 15 one-sided Monte Carlo p-values, with confirmation declared when `p_Fisher < 0.05`.

## Confirmatory result

Using 4,999 `k=1`-preserving Fourier surrogates per chromosome:

- Fisher statistic: `X = 104.456327692...`
- degrees of freedom: `30`
- nominal Fisher p-value: `3.63e-10`
- genome-level raw `p < 0.05`: `6/15`
- BH-significant: `6/15`
- Wilson 95% interval for the 6/15 proportion: approximately `19.8%–64.3%`
- resolution-stable genomes: `10/15`; `5/6` BH-significant genomes were stable.

Only one genome had zero surrogate exceedances and therefore attained the minimum reported Monte Carlo value `1/5000 = 0.0002`; no sub-grid p-value was imputed.

## Pre- versus post-unblinding materials

The unblinding boundary is explicit because it is central to the confirmatory design:

- `code/` and `protocols/` are **pre-unblinding** and were fixed before Confirmatory Cohort A sequence opening.
- `results/`, `figures/`, the confirmatory provenance outputs, and `docs/` are **post-unblinding outputs**.
- `tools/` contains **post-unblinding execution conveniences only**; these do not modify the frozen scientific analysis.

See `PROVENANCE_TIMELINE.md` and `FROZEN_CODE_SHA256SUMS`. A top-level `SHA256SUMS` records the release contents.

## Repository structure

- `code/` — exact pre-unblinding analysis scripts.
- `protocols/` — frozen preregistration and separately frozen tomography addendum.
- `manifests/` — accession lists, reserve queue, provenance hashes, and download metadata.
- `results/` — post-unblinding confirmatory result tables and summary statistics.
- `figures/` — post-unblinding manuscript figures.
- `docs/` — post-unblinding report and manuscript-support documents.
- `tools/` — additive post-unblinding execution wrapper for per-genome jobs.
- `release/` — hashes for key frozen bundles.

## Important confirmatory-status note

The tomography endpoint (`s4/s1`) was frozen in a separate addendum before sequence opening but was **not executed**, because the pre-existing tomography implementation was unavailable at unblinding. It was not reconstructed post hoc on the opened cohort. The addendum is retained here for transparency.

## Data acquisition

Raw chromosome FASTAs are not redistributed here. Download the exact nucleotide accessions listed in `manifests/CONFIRMATORY_COHORT_A.csv` from NCBI Nucleotide in FASTA format. The analysis validates accession/header and expected sequence length before computation.

## Reproduction

Tested environment for the archived run:

- Python 3.13.5
- NumPy 2.3.5
- pandas 2.2.3
- matplotlib 3.10.8

Install the dependencies in `requirements.txt`, place the 15 FASTAs in one directory, then run the canonical frozen runner:

```bash
python code/run_confirmatory_prevalence.py \
  --cohort manifests/CONFIRMATORY_COHORT_A.csv \
  --fasta-dir /path/to/confirmatory_fastas \
  --outdir reproduced_results \
  --null-reps 4999 \
  --seed 20260920
```

### Runtime benchmark

On the OpenAI Linux container used to package this release, two complete 4,999-surrogate single-genome benchmarks took **9.90 s** (`NC_017463.1`, 1.90 Mb) and **10.11 s** (`NZ_CP015506.1`, 5.46 Mb), with peak resident memory about **128 MB**. A serial 15-genome run would therefore be expected to take roughly **2.5–3 minutes on comparable hardware**, plus environment and I/O overhead. This is an illustrative benchmark, not a guaranteed runtime.

The analysis is not computationally heavyweight: on comparable hardware the full serial cohort is reproducible in roughly **under three minutes**. The earlier monolithic rerun failed because the interactive platform imposed a shorter single-command wall-clock limit, not because the analysis itself was expensive. For short-limit or CI environments, use the additive wrapper documented in `REPRODUCE.md`; it preserves the original per-genome seed schedule and leaves the frozen runner unchanged.

## Scientific scope

The confirmed result is a subset-limited computational observation about GC-skew phase architecture. It is **not** presented as a universal oriC-localization improvement, a new molecular mechanism, a thermodynamic identity, or a symbolic/metaphysical correspondence.

## License

Software and repository materials are released under the MIT License; see `LICENSE`.

## Citation

For the exact archived v1.0.0 release, cite:

Carretero R. *Symmetry-resolved GC-skew phase coherence: analysis and reproducibility package*. Version 1.0.0. Zenodo. 2026. DOI: **10.5281/zenodo.22909311**.

Repository-level concept DOI (all versions): **10.5281/zenodo.22909310**. The version-specific DOI above should be used when citing the exact archived release analyzed in the manuscript. See `CITATION.cff` and `ARCHIVAL_DOI.md`.
