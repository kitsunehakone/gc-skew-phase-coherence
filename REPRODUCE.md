# Reproduction checklist

## Canonical full-cohort path

1. Verify repository hashes with `sha256sum -c SHA256SUMS` where supported. At minimum, verify `FROZEN_CODE_SHA256SUMS`.
2. Obtain the 15 exact FASTA records listed in `manifests/CONFIRMATORY_COHORT_A.csv`.
3. Do not concatenate plasmids or other replicons.
4. Install Python dependencies from `requirements.txt`.
5. Run the canonical command shown in `README.md` with `--null-reps 4999 --seed 20260920`.
6. Compare the resulting per-genome values with `results/confirmatory_A_results.csv`.
7. Recompute Fisher's statistic directly from the stored per-genome `k1_preserved_phase_p` column:

```python
import numpy as np
import pandas as pd

df = pd.read_csv("results/confirmatory_A_results.csv")
p = df["k1_preserved_phase_p"].to_numpy(float)
print(-2 * np.log(p).sum())
# 104.456327692... expected
```

The Monte Carlo procedure is deterministic for the archived seed and frozen code. Minor differences in FFT/library behavior across platforms should be documented rather than silently adjusted.

## Short-wall-clock / CI path

`tools/run_confirmatory_chunked.py` is a **post-unblinding convenience wrapper**. It does not alter the statistic or null model. It invokes the frozen runner separately for each accession while preserving the exact original seed assignment

```text
seed_i = 20260920 + i * 100003
```

where `i` is the zero-based row index in the frozen cohort manifest.

Run genomes independently, for example:

```bash
python tools/run_confirmatory_chunked.py \
  --cohort manifests/CONFIRMATORY_COHORT_A.csv \
  --fasta-dir /path/to/confirmatory_fastas \
  --outdir chunked_reproduction \
  --null-reps 4999 \
  --seed 20260920 \
  --genome NC_017463.1
```

Repeat for each accession in the cohort manifest. Jobs may be distributed across machines because every genome has a deterministic independent seed assignment. After all 15 complete:

```bash
python tools/run_confirmatory_chunked.py \
  --cohort manifests/CONFIRMATORY_COHORT_A.csv \
  --outdir chunked_reproduction \
  --null-reps 4999 \
  --seed 20260920 \
  --merge
```

The merged files are written as `confirmatory_A_results_chunked.csv` and `confirmatory_A_summary_chunked.csv`. Compare their numerical contents with the archived `results/` tables.

## Runtime

Two full single-genome 4,999-surrogate benchmarks on the packaging environment took 9.90 s and 10.11 s, using approximately 128 MB peak resident memory. Expect runtime to vary with CPU, FFT library, filesystem and Python build.


## Chunked merge integrity

The wrapper is post-unblinding infrastructure only. It does not alter the statistic, null, cohort, or frozen code path. During merge it fails loudly if any accession is missing, if a result/provenance file is absent, if a result row lacks required fields or contains non-finite required values, if accessions disagree, or if a Monte Carlo p-value lies outside its allowed grid range.
