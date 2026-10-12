#!/usr/bin/env python3
from __future__ import annotations
import csv, json, math, os, wave
from pathlib import Path
import numpy as np

SOURCE_COMMIT = "4fa0cfcde7c919e1b68c99fdecd83ff9ac7334de"
B = 4999
BASE_SEED = 20261012
FMIN, FMAX = 70.0, 350.0
FRAME_SEC, HOP_SEC = 0.120, 0.060
ACF_MIN = 0.60
NH = 20
MIN_GOOD_FRAMES = 8
MIN_VALID_RECORDINGS = 25

FROZEN_FILES = [
"HC/002_i.wav","HC/004_i.wav","HC/006_i.wav","HC/016_i.wav","HC/024_i.wav",
"HC/026_i.wav","HC/028_i.wav","HC/042_i.wav","HC/049_i.wav","HC/053_i.wav",
"HC/061_i.wav","HC/063_i.wav","HC/065_i.wav","HC/077_i.wav","HC/081_i.wav",
"HC/085_i.wav","HC/086_i.wav","HC/089_i.wav","HC/097_i.wav","HC/099_i.wav",
"HC/107_i.wav","HC/109_i.wav","HC/111_i.wav","HC/113_i.wav","HC/115_i.wav",
"HC/117_i.wav","HC/119_i.wav","HC/121_i.wav","HC/123_i.wav","HC/125_i.wav",
"HC/127_i.wav","HC/129_i.wav","HC/131_i.wav",
]

def read_pcm_wav(path: Path):
    with wave.open(str(path), "rb") as w:
        nch = w.getnchannels()
        fs = w.getframerate()
        sw = w.getsampwidth()
        n = w.getnframes()
        comp = w.getcomptype()
        raw = w.readframes(n)
    if comp != "NONE":
        raise ValueError(f"compressed WAV: {comp}")
    if sw != 2:
        raise ValueError(f"expected 16-bit PCM, got sampwidth={sw}")
    x = np.frombuffer(raw, dtype="<i2").astype(np.float64) / 32768.0
    if nch != 1:
        raise ValueError(f"expected mono, got channels={nch}")
    return fs, x

def parabolic_peak(y, i):
    if i <= 0 or i >= len(y)-1:
        return float(i)
    a,b,c = float(y[i-1]), float(y[i]), float(y[i+1])
    den = a - 2*b + c
    if abs(den) < 1e-15:
        return float(i)
    delta = 0.5*(a-c)/den
    if abs(delta) > 1:
        return float(i)
    return float(i + delta)

def pitch_acf(frame, fs):
    x = np.asarray(frame, float)
    x = x - x.mean()
    e = float(np.dot(x,x))
    if e <= 1e-12:
        return None
    n = len(x)
    nfft = 1 << (2*n - 1).bit_length()
    X = np.fft.rfft(x, n=nfft)
    ac = np.fft.irfft(X*np.conj(X), n=nfft)[:n]
    # overlap correction makes the peak search less biased toward short lags
    lags = np.arange(n, dtype=float)
    ac = ac / np.maximum(1.0, n-lags)
    if ac[0] <= 0:
        return None
    acn = ac/ac[0]
    lo = max(1, int(math.floor(fs/FMAX)))
    hi = min(n-2, int(math.ceil(fs/FMIN)))
    if hi <= lo:
        return None
    i = lo + int(np.argmax(acn[lo:hi+1]))
    quality = float(acn[i])
    lag = parabolic_peak(acn, i)
    if lag <= 0:
        return None
    f0 = fs/lag
    if not (FMIN <= f0 <= FMAX):
        return None
    return f0, quality

def harmonic_coeffs(frame, fs, f0):
    x = np.asarray(frame, float) - np.mean(frame)
    h = np.hanning(len(x))
    t = np.arange(len(x), dtype=float)/fs
    norm = float(np.sum(h))
    coeff = np.empty(NH, dtype=np.complex128)
    for n in range(1, NH+1):
        coeff[n-1] = np.sum(h*x*np.exp(-2j*np.pi*n*f0*t))/norm
    return coeff

def recording_observed(path: Path):
    fs, x = read_pcm_wav(path)
    if fs != 44100:
        raise ValueError(f"sampling rate {fs} != 44100")
    if len(x) < fs:
        raise ValueError(f"duration {len(x)/fs:.3f}s < 1.0s")
    N = fs
    start = (len(x)-N)//2
    seg = np.array(x[start:start+N], float)
    seg -= seg.mean()
    L = int(round(FRAME_SEC*fs))
    H = int(round(HOP_SEC*fs))
    frames = []
    f0s = []
    qualities = []
    for s in range(0, len(seg)-L+1, H):
        fr = seg[s:s+L]
        pq = pitch_acf(fr, fs)
        if pq is None:
            continue
        f0, q = pq
        if q < ACF_MIN:
            continue
        c = harmonic_coeffs(fr, fs, f0)
        frames.append(c)
        f0s.append(f0)
        qualities.append(q)
    if len(frames) < MIN_GOOD_FRAMES:
        raise ValueError(f"accepted frames {len(frames)} < {MIN_GOOD_FRAMES}")
    C = np.vstack(frames)  # F x 20
    phase = np.angle(C)
    mag = np.abs(C)
    # harmonic centers n=2..19 map to zero-based indices j=1..18
    w = mag[:,0:18] * (mag[:,1:19]**2) * mag[:,2:20]
    kap = phase[:,2:20] + phase[:,0:18] - 2.0*phase[:,1:19]
    sw = float(w.sum())
    if sw <= 0:
        raise ValueError("zero closure weight")
    z = np.sum(w*np.exp(1j*kap))/sw
    neff = float(sw*sw/np.sum(w*w))
    return {
        "fs":fs, "accepted_frames":int(C.shape[0]),
        "median_f0":float(np.median(f0s)),
        "median_acf_quality":float(np.median(qualities)),
        "S1":float(np.real(z)), "M1_abs":float(abs(z)),
        "M1_angle_rad":float(np.angle(z)),
        "weight_neff":neff,
        "C":C, "w":w, "sw":sw,
    }

def null_p(obs, seed):
    C, w, sw = obs["C"], obs["w"], obs["sw"]
    F = C.shape[0]
    phi1 = np.angle(C[:,0])
    rng = np.random.default_rng(seed)
    null = np.empty(B, float)
    # Generate phases frame-by-frame/harmonic-by-harmonic under the frozen null.
    for b in range(B):
        ph = rng.uniform(0.0, 2*np.pi, size=(F,NH))
        ph[:,0] = phi1
        kap = ph[:,2:20] + ph[:,0:18] - 2.0*ph[:,1:19]
        null[b] = float(np.sum(w*np.cos(kap))/sw)
    p = float((1 + np.sum(null >= obs["S1"]))/(B+1))
    return p, {
        "null_mean":float(np.mean(null)),
        "null_sd":float(np.std(null,ddof=1)),
        "null_q95":float(np.quantile(null,0.95)),
        "null_q99":float(np.quantile(null,0.99)),
    }

def bh(pvals):
    p=np.asarray(pvals,float); o=np.argsort(p); r=p[o]; m=len(r)
    q=np.minimum.accumulate((r*m/np.arange(1,m+1))[::-1])[::-1]
    out=np.empty_like(q); out[o]=np.clip(q,0,1); return out

def chi2_sf_even(x,df):
    n=df//2; y=x/2
    return math.exp(-y)*sum(y**j/math.factorial(j) for j in range(n))

def main():
    root = Path(os.environ["SOUND_DATA_ROOT"])
    rows=[]
    for idx, rel in enumerate(sorted(FROZEN_FILES), start=1):
        path=root/rel
        row={"index":idx,"file":rel}
        try:
            obs=recording_observed(path)
            seed=BASE_SEED+100003*idx
            p,nd=null_p(obs,seed)
            row.update({
                "valid":True,
                "accepted_frames":obs["accepted_frames"],
                "median_f0":obs["median_f0"],
                "median_acf_quality":obs["median_acf_quality"],
                "S1":obs["S1"],
                "M1_abs":obs["M1_abs"],
                "M1_angle_rad":obs["M1_angle_rad"],
                "weight_neff":obs["weight_neff"],
                "p":p, **nd
            })
        except Exception as e:
            row.update({"valid":False,"reason":str(e)})
        rows.append(row)
        print("ROW",json.dumps(row,sort_keys=True),flush=True)

    valid=[r for r in rows if r["valid"]]
    pvals=[r["p"] for r in valid]
    if pvals:
        qvals=bh(pvals)
        for r,q in zip(valid,qvals): r["q_BH"]=float(q)
        X=float(-2*np.log(np.asarray(pvals)).sum())
        df=2*len(pvals)
        fp=chi2_sf_even(X,df)
    else:
        X=float("nan"); df=0; fp=float("nan")
    positive=(len(valid)>=MIN_VALID_RECORDINGS and fp<0.05)
    summary={
        "source_commit":SOURCE_COMMIT,
        "frozen_total":len(FROZEN_FILES),
        "valid_n":len(valid),
        "minimum_valid_required":MIN_VALID_RECORDINGS,
        "null_reps":B,
        "fisher_X":X,
        "fisher_df":df,
        "fisher_p":fp,
        "BLINDED_SOUND_I_REPLICATION_POSITIVE":bool(positive),
        "raw_p_lt_0_05_n":int(sum(r["p"]<0.05 for r in valid)),
        "BH_q_lt_0_05_n":int(sum(r.get("q_BH",1.0)<0.05 for r in valid)),
        "mean_S1":float(np.mean([r["S1"] for r in valid])) if valid else None,
        "median_S1":float(np.median([r["S1"] for r in valid])) if valid else None,
        "mean_M1_abs":float(np.mean([r["M1_abs"] for r in valid])) if valid else None,
        "median_F0":float(np.median([r["median_f0"] for r in valid])) if valid else None,
    }
    fields=sorted({k for r in rows for k in r.keys()})
    with open("sound_i_closure_results.csv","w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(rows)
    with open("sound_i_closure_summary.json","w") as f:
        json.dump(summary,f,indent=2,sort_keys=True)
    print("PRIMARY_SUMMARY",json.dumps(summary,sort_keys=True),flush=True)

if __name__=="__main__":
    main()
