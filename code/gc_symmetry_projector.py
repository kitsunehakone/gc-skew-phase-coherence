#!/usr/bin/env python3
"""
Symmetry-projector analysis for circular GC-skew signals.

Mathematical model
------------------
For an even-length real signal s[n], n in Z_M:
  H s[n] = s[n + M/2]
  R_a s[n] = s[a - n]       (mod M)

where a=0,...,M-1 indexes an unoriented reflection axis
phi_a = pi*a/M (mod pi).

The four orthogonal projectors are
  P_{r,h}^{(a)} = 1/4 (I + r R_a)(I + h H),
  r,h in {+1,-1}.

The target replication-sector statistic is
  J(a) = ||P_{-,-}^{(a)} s_c||^2 / ||s_c||^2,

where s_c is the centered GC-skew profile.

Efficient identity
------------------
Let u_- = (I-H)s_c / 2.  Then
  ||P_{-,-}^{(a)} s_c||^2
    = 1/2 ( ||u_-||^2 - <u_-, R_a u_-> ).

The quantity <u, R_a u> = sum_n u[n] u[a-n] is the circular
self-convolution (u*u)[a], so all M candidate axes are obtained
with a single FFT:
  c = ifft(fft(u) * fft(u)).real.

This makes the complete axis scan O(M log M).
"""

from __future__ import annotations
import argparse
import gzip
from pathlib import Path
import numpy as np


def read_fasta(path: str | Path) -> str:
    """Read one or more FASTA records and concatenate their sequence."""
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    chunks = []
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.startswith(">"):
                chunks.append(line.strip())
    seq = "".join(chunks).upper()
    if not seq:
        raise ValueError(f"No sequence found in {path}")
    return seq


def windowed_gc_skew(seq: str, bins: int = 512) -> np.ndarray:
    """
    Divide a circular sequence into `bins` approximately equal contiguous bins
    and compute (G-C)/(G+C) in each.  bins must be even.
    """
    if bins < 4 or bins % 2:
        raise ValueError("bins must be even and >= 4")
    n = len(seq)
    edges = np.linspace(0, n, bins + 1, dtype=int)
    out = np.empty(bins, dtype=float)
    for i in range(bins):
        w = seq[edges[i]:edges[i+1]]
        g = w.count("G")
        c = w.count("C")
        d = g + c
        out[i] = (g - c) / d if d else 0.0
    return out


def _self_convolution_real(x: np.ndarray) -> np.ndarray:
    """Circular convolution x*x for real x."""
    f = np.fft.fft(x)
    return np.fft.ifft(f * f).real


def sector_scan(signal: np.ndarray):
    """
    Return exact four-sector energies for every reflection axis a.

    Keys are '++', '-+', '+-', '--' where first sign is R parity,
    second sign is H parity.  Each value has shape (M,).
    """
    s = np.asarray(signal, dtype=float)
    M = len(s)
    if M < 4 or M % 2:
        raise ValueError("signal length must be even and >= 4")
    s = s - s.mean()

    Hs = np.roll(s, -M // 2)
    hp = 0.5 * (s + Hs)
    hm = 0.5 * (s - Hs)

    ep = float(np.dot(hp, hp))
    em = float(np.dot(hm, hm))
    cp = _self_convolution_real(hp)
    cm = _self_convolution_real(hm)

    # On each H eigenspace, (I +/- R)/2 splits energy orthogonally.
    e_pp = 0.5 * (ep + cp)  # R+, H+
    e_mp = 0.5 * (ep - cp)  # R-, H+
    e_pm = 0.5 * (em + cm)  # R+, H-
    e_mm = 0.5 * (em - cm)  # R-, H-

    # Remove tiny roundoff negatives only.
    tol = 1e-10 * max(1.0, float(np.dot(s, s)))
    arrays = {}
    for key, x in {"++": e_pp, "-+": e_mp, "+-": e_pm, "--": e_mm}.items():
        x = np.where((x < 0) & (x > -tol), 0.0, x)
        arrays[key] = x

    return s, arrays


def scan_J(signal: np.ndarray):
    """Return J(a), sector-energy fractions, and entropy for every axis."""
    s, energies = sector_scan(signal)
    total = float(np.dot(s, s))
    if total <= 0:
        raise ValueError("centered signal is zero")

    fractions = {k: v / total for k, v in energies.items()}
    J = fractions["--"]

    # Entropy of four (R,H) sector fractions at each axis.
    P = np.vstack([fractions[k] for k in ("++", "-+", "+-", "--")]).T
    P = np.clip(P, 0.0, 1.0)
    row_sums = P.sum(axis=1)
    P = P / row_sums[:, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        logs = np.where(P > 0, np.log2(P), 0.0)
    entropy = -(P * logs).sum(axis=1)

    return J, fractions, entropy


def best_projector_axis(signal: np.ndarray):
    """
    Return best discrete axis.
    a in {0,...,M-1}; phi = pi*a/M in [0,pi).
    """
    J, fractions, entropy = scan_J(signal)
    a = int(np.argmax(J))
    M = len(J)
    return {
        "axis_index": a,
        "phi_rad": np.pi * a / M,
        "phi_fraction_of_circle": a / (2.0 * M),
        "J": float(J[a]),
        "entropy_bits": float(entropy[a]),
        "fractions": {k: float(v[a]) for k, v in fractions.items()},
    }


def first_harmonic_axis(signal: np.ndarray):
    """
    First-harmonic baseline.

    If s(theta) ~ A sin(theta - phi), then arg(F_1) = -phi-pi/2.
    Axis is unoriented, so phi is reduced modulo pi.
    """
    s = np.asarray(signal, dtype=float)
    s = s - s.mean()
    M = len(s)
    theta = 2 * np.pi * np.arange(M) / M
    F1 = np.dot(s, np.exp(-1j * theta))
    if abs(F1) == 0:
        return {"phi_rad": np.nan, "power_fraction": 0.0}

    phi = (-np.angle(F1) - np.pi / 2.0) % np.pi

    # Energy fraction in k=±1 pair, using an orthonormal FFT.
    F = np.fft.fft(s, norm="ortho")
    total = float(np.sum(np.abs(F) ** 2))
    p1 = float((abs(F[1]) ** 2 + abs(F[-1]) ** 2) / total) if total else 0.0
    return {"phi_rad": float(phi), "power_fraction": p1}


def cumulative_skew_axis(signal: np.ndarray):
    """
    Simple cumulative-skew baseline: return the axis through the cumulative
    minimum.  Because the axis is unoriented, phi and phi+pi are equivalent.
    """
    s = np.asarray(signal, dtype=float)
    s = s - s.mean()
    c = np.cumsum(s)
    i = int(np.argmin(c))
    M = len(s)
    phi = (2 * np.pi * i / M) % np.pi
    return {"phi_rad": float(phi), "index": i}


def axis_error(phi_hat: float, phi_true: float) -> float:
    """Unoriented angular-axis error in radians, range [0, pi/2]."""
    if not np.isfinite(phi_hat):
        return np.nan
    d = (phi_hat - phi_true) % np.pi
    return float(min(d, np.pi - d))


def analyze_fasta(path: str | Path, bins: int = 512):
    seq = read_fasta(path)
    skew = windowed_gc_skew(seq, bins=bins)
    proj = best_projector_axis(skew)
    fh = first_harmonic_axis(skew)
    cum = cumulative_skew_axis(skew)
    return {
        "sequence_length": len(seq),
        "bins": bins,
        "projector": proj,
        "first_harmonic": fh,
        "cumulative": cum,
        "skew": skew,
    }


def _self_test():
    # Exact odd-sine signal about a chosen axis.
    M = 512
    phi0 = 0.731
    th = 2 * np.pi * np.arange(M) / M
    s = (
        0.20 * np.sin(th - phi0)
        + 1.00 * np.sin(3 * (th - phi0))
        + 0.65 * np.sin(5 * (th - phi0))
        - 0.35 * np.sin(9 * (th - phi0))
    )

    result = best_projector_axis(s)
    err = axis_error(result["phi_rad"], phi0)
    grid = np.pi / M
    assert result["J"] > 0.999, result
    assert err <= grid + 1e-12, (err, grid, result)

    # Parseval/projector identity at every axis.
    sc, E = sector_scan(s)
    total = float(np.dot(sc, sc))
    esum = E["++"] + E["-+"] + E["+-"] + E["--"]
    assert np.max(np.abs(esum - total)) < 1e-8 * total

    print("SELF-TEST PASS")
    print(f"true axis phi={phi0:.6f}")
    print(f"estimated phi={result['phi_rad']:.6f}")
    print(f"J={result['J']:.12f}")
    print(f"entropy={result['entropy_bits']:.6g} bits")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fasta", nargs="?", help="FASTA or FASTA.gz genome")
    ap.add_argument("--bins", type=int, default=512)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()

    if args.self_test:
        _self_test()
        return

    if not args.fasta:
        ap.error("provide FASTA path or use --self-test")

    r = analyze_fasta(args.fasta, bins=args.bins)
    p = r["projector"]
    fh = r["first_harmonic"]
    cu = r["cumulative"]
    print(f"length_bp={r['sequence_length']}")
    print(f"bins={r['bins']}")
    print(f"projector_axis_phi_rad={p['phi_rad']:.12g}")
    print(f"projector_axis_fraction_circle={p['phi_fraction_of_circle']:.12g}")
    print(f"J={p['J']:.12g}")
    print(f"symmetry_entropy_bits={p['entropy_bits']:.12g}")
    for k, v in p["fractions"].items():
        print(f"sector_{k}_fraction={v:.12g}")
    print(f"first_harmonic_phi_rad={fh['phi_rad']:.12g}")
    print(f"first_harmonic_power_fraction={fh['power_fraction']:.12g}")
    print(f"cumulative_phi_rad={cu['phi_rad']:.12g}")


if __name__ == "__main__":
    main()
