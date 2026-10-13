import json, math, os, re, sys
from pathlib import Path

import numpy as np
from scipy.io import wavfile
from scipy import signal, optimize, special

SOURCE_COMMIT = "4fa0cfcde7c919e1b68c99fdecd83ff9ac7334de"
CAL_IDS = ["008","021","024","027","031","039","046","052","058","064","072","078","084","094","098","102"]
VAL_IDS = ["020","022","025","028","032","042","048","055","062","068","076","080","092","096","100"]

FRAME_SEC = 0.120
HOP_SEC = 0.060
CENTER_SEC = 1.0
F0_MIN = 70.0
F0_MAX = 350.0
ACF_MIN = 0.60
GRID_N = 512

def wrap(x):
    return np.angle(np.exp(1j*x))

def read_mono(path):
    fs, x = wavfile.read(path)
    if x.ndim > 1:
        x = x.mean(axis=1)
    if np.issubdtype(x.dtype, np.integer):
        info = np.iinfo(x.dtype)
        scale = max(abs(info.min), abs(info.max))
        x = x.astype(np.float64) / float(scale)
    else:
        x = x.astype(np.float64)
    return int(fs), x

def f0_from_acf(frame, fs):
    x = np.asarray(frame, float)
    x = x - np.mean(x)
    e = float(np.dot(x, x))
    if e <= 1e-14:
        return None, 0.0
    ac = signal.fftconvolve(x, x[::-1], mode="full")[len(x)-1:]
    ac = ac / ac[0]
    lag_lo = max(1, int(math.floor(fs/F0_MAX)))
    lag_hi = min(len(ac)-1, int(math.ceil(fs/F0_MIN)))
    if lag_hi <= lag_lo:
        return None, 0.0
    seg = ac[lag_lo:lag_hi+1]
    rel = int(np.argmax(seg))
    lag = lag_lo + rel
    peak = float(ac[lag])
    if peak < ACF_MIN:
        return None, peak
    # Parabolic interpolation of autocorrelation peak where possible.
    lag_f = float(lag)
    if 1 <= lag < len(ac)-1:
        y0,y1,y2 = ac[lag-1],ac[lag],ac[lag+1]
        den = y0 - 2*y1 + y2
        if abs(den) > 1e-12:
            delta = 0.5*(y0-y2)/den
            if abs(delta) <= 1:
                lag_f += float(delta)
    return fs/lag_f, peak

def frame_closures(frame, fs):
    f0, peak = f0_from_acf(frame, fs)
    if f0 is None:
        return None
    if 4*f0 >= 0.9*(fs/2):
        return None
    n = len(frame)
    w = np.hanning(n)
    xw = (frame - np.mean(frame))*w
    t = np.arange(n)/fs
    cs = []
    for h in range(1,5):
        c = np.sum(xw*np.exp(-2j*np.pi*h*f0*t))
        cs.append(c)
    if min(abs(c) for c in cs) <= 1e-12:
        return None
    ph = np.angle(cs)
    k1 = wrap(ph[0] - 2*ph[1] + ph[2])
    k2 = wrap(ph[1] - 2*ph[2] + ph[3])
    chi = wrap(k1+k2)
    return float(k1), float(k2), float(chi), float(f0), float(peak)

def process_record(path):
    fs,x = read_mono(path)
    n_center = int(round(CENTER_SEC*fs))
    if len(x) < n_center:
        return {"valid": False, "reason": "short_recording", "frames": []}
    mid = len(x)//2
    start = mid - n_center//2
    seg = x[start:start+n_center]
    L = int(round(FRAME_SEC*fs))
    H = int(round(HOP_SEC*fs))
    rows = []
    for s in range(0, len(seg)-L+1, H):
        out = frame_closures(seg[s:s+L], fs)
        if out is not None:
            k1,k2,chi,f0,peak = out
            rows.append({"k1":k1,"k2":k2,"chi":chi,"f0":f0,"acf_peak":peak})
    return {"valid": len(rows)>=8, "reason": None if len(rows)>=8 else "too_few_frames", "frames": rows, "fs":fs}

# Deterministic torus grid.
grid = np.linspace(-np.pi, np.pi, GRID_N, endpoint=False)
g1,g2 = np.meshgrid(grid,grid,indexing="ij")
ct1 = np.cos(g1)
ct2 = np.cos(g2)
ct3 = np.cos(g1+g2)
Tgrid = ct1+ct2+2*ct3
log_area_per = 2*math.log(2*math.pi) - 2*math.log(GRID_N)

def logZ_one(j):
    return float(special.logsumexp(j*Tgrid) + log_area_per)

def moments_one(j):
    a = j*Tgrid
    m = np.max(a)
    w = np.exp(a-m)
    den = w.sum()
    return np.array([(w*ct1).sum()/den,(w*ct2).sum()/den,(w*ct3).sum()/den], float)

def fit_one(k1,k2):
    T = np.cos(k1)+np.cos(k2)+2*np.cos(k1+k2)
    def nll(j):
        return -(j*T.sum() - len(T)*logZ_one(float(j)))
    opt = optimize.minimize_scalar(nll,bounds=(0.0,20.0),method="bounded",options={"xatol":1e-10})
    return float(opt.x), float(opt.fun)

def logZ_three(par):
    j1,j2,j3 = par
    a = j1*ct1+j2*ct2+j3*ct3
    return float(special.logsumexp(a)+log_area_per)

def fit_three(k1,k2):
    c1=np.cos(k1); c2=np.cos(k2); c3=np.cos(k1+k2)
    def nll(par):
        j1,j2,j3=par
        return -(j1*c1.sum()+j2*c2.sum()+j3*c3.sum()-len(k1)*logZ_three(par))
    opt=optimize.minimize(nll,x0=np.array([0.1,0.1,0.2]),method="L-BFGS-B",bounds=[(-20,20)]*3)
    return [float(x) for x in opt.x], float(opt.fun)

def circ_stats(x):
    z=np.mean(np.exp(1j*np.asarray(x)))
    return {"mean_angle":float(np.angle(z)),"resultant":float(abs(z))}

root=Path(os.environ.get("MINSK_SOURCE","source"))
all_ids=CAL_IDS+VAL_IDS
records={}
for sid in all_ids:
    p=root/"ALS"/f"{sid}_a.wav"
    if not p.exists():
        records[sid]={"valid":False,"reason":"missing","frames":[]}
    else:
        records[sid]=process_record(p)

def collect(ids):
    rows=[]
    valid_ids=[]
    for sid in ids:
        r=records[sid]
        if r.get("valid"):
            valid_ids.append(sid)
            for f in r["frames"]:
                rows.append((sid,f["k1"],f["k2"],f["chi"],f["f0"],f["acf_peak"]))
    return valid_ids,rows

cal_valid,cal_rows=collect(CAL_IDS)
val_valid,val_rows=collect(VAL_IDS)

if len(cal_valid)==0 or len(val_valid)==0:
    raise RuntimeError("No valid calibration or validation recordings")

cal=np.array([[r[1],r[2],r[3]] for r in cal_rows],float)
val=np.array([[r[1],r[2],r[3]] for r in val_rows],float)

j,nll=fit_one(cal[:,0],cal[:,1])
pred=moments_one(j)
emp=np.array([np.cos(val[:,0]).mean(),np.cos(val[:,1]).mean(),np.cos(val[:,2]).mean()])
disc=float(np.linalg.norm(emp-pred))
logz=logZ_one(j)
Tval=np.cos(val[:,0])+np.cos(val[:,1])+2*np.cos(val[:,2])
lval=float(np.mean(j*Tval-logz))
uniform=-2*math.log(2*math.pi)

# participant-clustered moment means
cluster=[]
for sid in val_valid:
    a=np.array([[f["k1"],f["k2"],f["chi"]] for f in records[sid]["frames"]])
    cluster.append([np.cos(a[:,0]).mean(),np.cos(a[:,1]).mean(),np.cos(a[:,2]).mean()])
cluster=np.array(cluster,float)
cluster_mean=cluster.mean(axis=0)

# unconstrained diagnostic
j3,nll3=fit_three(cal[:,0],cal[:,1])
logz3=logZ_three(j3)
lval3=float(np.mean(j3[0]*np.cos(val[:,0])+j3[1]*np.cos(val[:,1])+j3[2]*np.cos(val[:,2])-logz3))

# stiffness/variance diagnostic using wrapped +/- coordinates on validation
kp=wrap((val[:,0]+val[:,1])/math.sqrt(2))
km=wrap((val[:,0]-val[:,1])/math.sqrt(2))
varp=float(np.mean(kp**2))
varm=float(np.mean(km**2))
var_ratio=float(varm/varp) if varp>0 else None

success = bool(
    j > 1e-8
    and lval > uniform
    and emp[2] > emp[0]
    and emp[2] > emp[1]
    and disc <= 0.15
)

result={
    "source_commit":SOURCE_COMMIT,
    "calibration_ids":CAL_IDS,
    "validation_ids":VAL_IDS,
    "valid_calibration_recordings":cal_valid,
    "valid_validation_recordings":val_valid,
    "calibration_frames":int(len(cal)),
    "validation_frames":int(len(val)),
    "invalid_recordings":{sid:records[sid].get("reason") for sid in all_ids if not records[sid].get("valid")},
    "fitted_j":j,
    "predicted_moments":{"M1":float(pred[0]),"M2":float(pred[1]),"Mchi":float(pred[2])},
    "validation_moments":{"M1":float(emp[0]),"M2":float(emp[1]),"Mchi":float(emp[2])},
    "validation_clustered_moments":{"M1":float(cluster_mean[0]),"M2":float(cluster_mean[1]),"Mchi":float(cluster_mean[2])},
    "moment_discrepancy":disc,
    "validation_mean_log_density":lval,
    "uniform_mean_log_density":uniform,
    "validation_log_density_gain":lval-uniform,
    "primary_success":success,
    "unconstrained_fit":{"j1":j3[0],"j2":j3[1],"j3":j3[2],"validation_mean_log_density":lval3,"gain_vs_k4":lval3-lval},
    "circular_validation":{"kappa1":circ_stats(val[:,0]),"kappa2":circ_stats(val[:,1]),"chi":circ_stats(val[:,2])},
    "small_fluctuation_diagnostic":{"var_kappa_plus":varp,"var_kappa_minus":varm,"var_minus_over_plus":var_ratio,"K4_target_if_small_fluctuation":5.0},
    "per_recording":{sid:{"valid":records[sid]["valid"],"n_frames":len(records[sid].get("frames",[])),
                         "median_f0":float(np.median([f["f0"] for f in records[sid]["frames"]])) if records[sid].get("frames") else None}
                     for sid in all_ids}
}

Path("results").mkdir(exist_ok=True)
Path("results/k4_harmonic_als_a_results.json").write_text(json.dumps(result,indent=2))

lines=[]
lines.append("K4 harmonic coherence prospective ALS /a/ test")
lines.append(f"source commit: {SOURCE_COMMIT}")
lines.append(f"valid calibration recordings: {len(cal_valid)}/{len(CAL_IDS)}; frames={len(cal)}")
lines.append(f"valid validation recordings: {len(val_valid)}/{len(VAL_IDS)}; frames={len(val)}")
lines.append(f"fitted j: {j:.12g}")
lines.append(f"predicted moments: M1={pred[0]:.8f} M2={pred[1]:.8f} Mchi={pred[2]:.8f}")
lines.append(f"validation moments: M1={emp[0]:.8f} M2={emp[1]:.8f} Mchi={emp[2]:.8f}")
lines.append(f"clustered validation: M1={cluster_mean[0]:.8f} M2={cluster_mean[1]:.8f} Mchi={cluster_mean[2]:.8f}")
lines.append(f"moment discrepancy: {disc:.8f}")
lines.append(f"held-out mean log density: {lval:.8f}")
lines.append(f"uniform mean log density: {uniform:.8f}")
lines.append(f"log-density gain: {lval-uniform:.8f}")
lines.append(f"PRIMARY SUCCESS: {success}")
lines.append(f"unconstrained calibration fit j1:j2:j3 = {j3[0]:.8f}:{j3[1]:.8f}:{j3[2]:.8f}")
lines.append(f"unconstrained held-out mean log density: {lval3:.8f}; gain vs K4={lval3-lval:.8f}")
lines.append(f"small-fluctuation var_minus/var_plus={var_ratio:.8f} (K4 target 5 if regime applicable)")
Path("results/k4_harmonic_als_a_results.txt").write_text("\n".join(lines)+"\n")
print("\n".join(lines))
