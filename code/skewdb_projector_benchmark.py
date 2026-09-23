#!/usr/bin/env python3
"""
Batch benchmark pipeline for the GC-skew symmetry projector.

Inputs
------
1. SkewDB skplot.csv (or any subset of it), whose relevant fields are:
   name, relpos, abspos, gcskew, gccount
   where gcskew and gccount are cumulative quantities.
2. Optional SkewDB gcskewdb.csv summary for metadata / genome sizes.
3. Optional truth CSV (e.g. DoriC export) containing an accession and oriC
   start/end coordinates.

Outputs
-------
One row per sequence with:
  J_max
  symmetry entropy
  four (R,H) sector fractions at the best axis
  first-harmonic power fraction
  half-turn-odd energy fraction
  higher-odd energy fraction
  odd-axis coherence
  J gain over k=1
  projector / first-harmonic / cumulative axis estimates
  optional errors to oriC ground truth
  optional phase-randomization p-values

The key null preserves the complete Fourier magnitude spectrum but randomizes
phases.  A stricter version preserves k=1 exactly and randomizes only the
higher Fourier phases.  Thus a significant result tests multi-harmonic phase
coherence beyond first-harmonic information.

This script expects gc_symmetry_projector.py to be in the same directory.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
import numpy as np
import pandas as pd

from gc_symmetry_projector import (
    best_projector_axis,
    first_harmonic_axis,
    cumulative_skew_axis,
    scan_J,
    axis_error,
)


def local_skew_from_cumulative(group: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """
    Reconstruct local normalized GC skew from SkewDB cumulative quantities.

    If C_j is cumulative (G-C) and Q_j cumulative (G+C), then on interval j:
        local_skew_j = (C_j-C_{j-1}) / (Q_j-Q_{j-1}).

    Returns
    -------
    centers : relative positions in [0,1)
    skew    : local normalized GC skew
    """
    g = group.sort_values("abspos").copy()
    gc = pd.to_numeric(g["gcskew"], errors="coerce").to_numpy(float)
    gcc = pd.to_numeric(g["gccount"], errors="coerce").to_numpy(float)

    if "relpos" in g.columns:
        endpoints = pd.to_numeric(g["relpos"], errors="coerce").to_numpy(float)
        # Some exports may use percentages rather than fractions.
        finite = endpoints[np.isfinite(endpoints)]
        if len(finite) and np.nanmax(finite) > 1.5:
            endpoints = endpoints / np.nanmax(endpoints)
    else:
        apos = pd.to_numeric(g["abspos"], errors="coerce").to_numpy(float)
        endpoints = apos / np.nanmax(apos)

    dgc = np.diff(np.r_[0.0, gc])
    dgcc = np.diff(np.r_[0.0, gcc])
    with np.errstate(divide="ignore", invalid="ignore"):
        skew = np.where(dgcc > 0, dgc / dgcc, np.nan)

    starts = np.r_[0.0, endpoints[:-1]]
    centers = 0.5 * (starts + endpoints)

    mask = np.isfinite(centers) & np.isfinite(skew)
    centers = centers[mask] % 1.0
    skew = skew[mask]

    order = np.argsort(centers)
    return centers[order], skew[order]


def periodic_resample(x: np.ndarray, y: np.ndarray, bins: int = 512) -> np.ndarray:
    """Periodically interpolate an irregular circular profile onto an even grid."""
    if bins < 4 or bins % 2:
        raise ValueError("bins must be even and >=4")
    x = np.asarray(x, float) % 1.0
    y = np.asarray(y, float)
    if len(x) < 4:
        raise ValueError("too few valid local-skew bins")

    order = np.argsort(x)
    x, y = x[order], y[order]

    # Merge duplicate positions, if any.
    ux, inv = np.unique(x, return_inverse=True)
    if len(ux) != len(x):
        sums = np.zeros(len(ux))
        counts = np.zeros(len(ux))
        np.add.at(sums, inv, y)
        np.add.at(counts, inv, 1)
        x, y = ux, sums / counts

    xe = np.r_[x[-1] - 1.0, x, x[0] + 1.0]
    ye = np.r_[y[-1], y, y[0]]
    target = (np.arange(bins) + 0.5) / bins
    return np.interp(target, xe, ye)


def spectral_energy_stats(signal: np.ndarray) -> dict:
    """Fourier quantities using an orthonormal FFT."""
    s = np.asarray(signal, float)
    s = s - s.mean()
    M = len(s)
    F = np.fft.fft(s, norm="ortho")
    e = np.abs(F) ** 2
    total = float(e.sum())
    if total <= 0:
        return dict(
            first_harmonic_power_fraction=0.0,
            hminus_energy_fraction=0.0,
            higher_odd_energy_fraction=0.0,
        )

    p1 = float((e[1] + e[-1]) / total)

    # For even M, H=T^(M/2), so H=-1 exactly on odd k.
    odd_mask = (np.arange(M) % 2) == 1
    hminus = float(e[odd_mask].sum() / total)
    higher_odd = max(0.0, hminus - p1)
    return dict(
        first_harmonic_power_fraction=p1,
        hminus_energy_fraction=hminus,
        higher_odd_energy_fraction=higher_odd,
    )


def randomize_fourier_phases(
    signal: np.ndarray,
    rng: np.random.Generator,
    preserve_k1: bool = False,
) -> np.ndarray:
    """
    Preserve the exact Fourier magnitude spectrum of a real signal while
    randomizing phases.  If preserve_k1=True, the k=1 conjugate pair is kept
    exactly, isolating information contained in higher-harmonic phase coherence.
    """
    s = np.asarray(signal, float)
    s = s - s.mean()
    M = len(s)
    if M % 2:
        raise ValueError("signal length must be even")

    F = np.fft.rfft(s)
    G = np.zeros_like(F, dtype=complex)
    G[0] = F[0]  # centered signal: normally zero

    ny = M // 2
    for k in range(1, ny):
        if preserve_k1 and k == 1:
            G[k] = F[k]
        else:
            mag = abs(F[k])
            phase = rng.uniform(0.0, 2.0 * np.pi)
            G[k] = mag * np.exp(1j * phase)

    # Nyquist coefficient of a real signal must remain real.
    G[ny] = F[ny].real
    return np.fft.irfft(G, n=M)


def phase_coherence_test(
    signal: np.ndarray,
    n_null: int = 499,
    seed: int = 20260919,
    preserve_k1: bool = True,
) -> dict:
    """
    Monte Carlo test of observed J_max against phase-randomized surrogates.

    Null preserves:
      * the entire Fourier magnitude spectrum;
      * optionally the complete k=1 coefficient (magnitude + phase).

    p = (1 + #{J_null >= J_obs}) / (n_null + 1)
    """
    obs = best_projector_axis(signal)["J"]
    rng = np.random.default_rng(seed)
    null = np.empty(n_null, float)
    for i in range(n_null):
        sur = randomize_fourier_phases(signal, rng, preserve_k1=preserve_k1)
        null[i] = best_projector_axis(sur)["J"]

    p = (1.0 + np.sum(null >= obs)) / (n_null + 1.0)
    return {
        "J_observed": float(obs),
        "null_mean_J": float(np.mean(null)),
        "null_q95_J": float(np.quantile(null, 0.95)),
        "null_q99_J": float(np.quantile(null, 0.99)),
        "phase_coherence_p": float(p),
    }


def circ_axis_fraction(phi_rad: float) -> float:
    """Convert unoriented axis phi in [0,pi) to sequence fraction in [0,0.5)."""
    if not np.isfinite(phi_rad):
        return np.nan
    return float((phi_rad % np.pi) / (2.0 * np.pi))


def ori_axis_error_fraction(pred_axis_frac: float, ori_frac: float) -> float:
    """
    Error between an unoriented axis and an oriC point, as fraction of full
    chromosome circumference, in [0,0.25].
    """
    a = pred_axis_frac % 0.5
    b = ori_frac % 0.5
    d = abs(a - b) % 0.5
    return float(min(d, 0.5 - d))


def analyze_signal(signal: np.ndarray, n_null: int = 0, seed: int = 20260919) -> dict:
    proj = best_projector_axis(signal)
    f1 = first_harmonic_axis(signal)
    cum = cumulative_skew_axis(signal)
    spec = spectral_energy_stats(signal)

    hminus = spec["hminus_energy_fraction"]
    coherence = proj["J"] / hminus if hminus > 0 else np.nan

    row = {
        "J_max": proj["J"],
        "symmetry_entropy_bits": proj["entropy_bits"],
        "E_pp_fraction": proj["fractions"]["++"],
        "E_mp_fraction": proj["fractions"]["-+"],
        "E_pm_fraction": proj["fractions"]["+-"],
        "E_mm_fraction": proj["fractions"]["--"],
        "first_harmonic_power_fraction": spec["first_harmonic_power_fraction"],
        "hminus_energy_fraction": hminus,
        "higher_odd_energy_fraction": spec["higher_odd_energy_fraction"],
        "odd_axis_coherence": coherence,
        "J_gain_over_k1": proj["J"] - spec["first_harmonic_power_fraction"],
        "projector_axis_fraction": proj["phi_fraction_of_circle"],
        "first_harmonic_axis_fraction": circ_axis_fraction(f1["phi_rad"]),
        "cumulative_axis_fraction": circ_axis_fraction(cum["phi_rad"]),
    }

    if n_null > 0:
        # Primary test: preserve k=1 exactly and randomize higher phases.
        nk1 = phase_coherence_test(
            signal, n_null=n_null, seed=seed, preserve_k1=True
        )
        row.update({
            "k1_preserved_null_mean_J": nk1["null_mean_J"],
            "k1_preserved_null_q95_J": nk1["null_q95_J"],
            "k1_preserved_phase_p": nk1["phase_coherence_p"],
        })

        # Secondary test: randomize every nontrivial Fourier phase.
        nall = phase_coherence_test(
            signal, n_null=n_null, seed=seed + 1, preserve_k1=False
        )
        row.update({
            "allphase_null_mean_J": nall["null_mean_J"],
            "allphase_null_q95_J": nall["null_q95_J"],
            "allphase_phase_p": nall["phase_coherence_p"],
        })

    return row


def add_truth_errors(
    row: dict,
    ori_start: float,
    ori_end: float,
    genome_size: float,
) -> dict:
    """Add axis errors to the midpoint of an annotated oriC interval."""
    if not (np.isfinite(ori_start) and np.isfinite(ori_end) and genome_size > 0):
        return row
    mid = 0.5 * (ori_start + ori_end)
    ori_frac = (mid / genome_size) % 1.0
    row["ori_mid_fraction"] = ori_frac
    for method in ("projector", "first_harmonic", "cumulative"):
        a = row[f"{method}_axis_fraction"]
        ef = ori_axis_error_fraction(a, ori_frac)
        row[f"{method}_axis_error_fraction"] = ef
        row[f"{method}_axis_error_degrees"] = ef * 360.0
        row[f"{method}_axis_error_bp"] = ef * genome_size
    return row


def load_truth(
    path: str | None,
    id_col: str,
    start_col: str,
    end_col: str,
    size_col: str | None = None,
) -> dict:
    """Load a generic oriC truth file keyed by accession."""
    if not path:
        return {}
    df = pd.read_csv(path)
    out = {}
    for _, r in df.iterrows():
        key = str(r[id_col])
        out[key] = {
            "start": float(r[start_col]),
            "end": float(r[end_col]),
            "size": float(r[size_col]) if size_col and size_col in df.columns else np.nan,
        }
    return out


def iter_skplot_groups(path: str, names: set[str] | None = None, chunksize: int = 1_000_000):
    """
    Stream a name-sorted skplot.csv without loading the multi-GB database.

    Assumption: rows for each `name` are contiguous, as in the SkewDB export.
    If using a custom unsorted CSV, sort by name, abspos first.
    """
    usecols = ["name", "relpos", "abspos", "gcskew", "gccount"]
    carry = None

    for chunk in pd.read_csv(path, usecols=usecols, chunksize=chunksize):
        if carry is not None:
            chunk = pd.concat([carry, chunk], ignore_index=True)
            carry = None

        if chunk.empty:
            continue

        last_name = chunk["name"].iloc[-1]
        carry = chunk[chunk["name"] == last_name].copy()
        done = chunk[chunk["name"] != last_name]

        for name, g in done.groupby("name", sort=False):
            if names is None or str(name) in names:
                yield str(name), g

    if carry is not None and not carry.empty:
        for name, g in carry.groupby("name", sort=False):
            if names is None or str(name) in names:
                yield str(name), g


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skplot", required=True, help="SkewDB skplot.csv or subset")
    ap.add_argument("--output", required=True)
    ap.add_argument("--accessions", help="text file, one sequence accession per line")
    ap.add_argument("--bins", type=int, default=512)
    ap.add_argument("--null-reps", type=int, default=0,
                    help="phase-null replicates per sequence; use 499 or 999 for final analysis")
    ap.add_argument("--seed", type=int, default=20260919)

    # Generic truth table mapping, compatible with a manually normalized DoriC export.
    ap.add_argument("--truth")
    ap.add_argument("--truth-id-col", default="refseq")
    ap.add_argument("--truth-start-col", default="oric_start")
    ap.add_argument("--truth-end-col", default="oric_end")
    ap.add_argument("--truth-size-col", default="genome_size")

    args = ap.parse_args()

    names = None
    if args.accessions:
        names = {
            x.strip() for x in Path(args.accessions).read_text().splitlines()
            if x.strip()
        }

    truth = load_truth(
        args.truth,
        args.truth_id_col,
        args.truth_start_col,
        args.truth_end_col,
        args.truth_size_col,
    )

    rows = []
    for i, (name, g) in enumerate(iter_skplot_groups(args.skplot, names=names)):
        try:
            x, y = local_skew_from_cumulative(g)
            signal = periodic_resample(x, y, bins=args.bins)
            row = {"name": name, "raw_points": len(g), "analysis_bins": args.bins}
            row.update(analyze_signal(
                signal,
                n_null=args.null_reps,
                seed=args.seed + i * 100003,
            ))

            if name in truth:
                t = truth[name]
                genome_size = t["size"]
                if not np.isfinite(genome_size):
                    genome_size = float(pd.to_numeric(g["abspos"], errors="coerce").max())
                row = add_truth_errors(row, t["start"], t["end"], genome_size)

            rows.append(row)
        except Exception as e:
            rows.append({"name": name, "error": repr(e)})

    out = pd.DataFrame(rows)
    out.to_csv(args.output, index=False)
    print(f"Wrote {len(out)} rows to {args.output}")


if __name__ == "__main__":
    main()
