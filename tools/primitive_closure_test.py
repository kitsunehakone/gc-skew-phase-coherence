#!/usr/bin/env python3
from __future__ import annotations
import csv, json, math, sys, time, urllib.parse, urllib.request
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
from gc_symmetry_projector import windowed_gc_skew, best_projector_axis

N_NULL = 4999
BASE_SEED = 20260920
BINS_LIST = (256, 512, 1024)

def fetch_fasta(acc: str, tries: int = 4) -> str:
    params = urllib.parse.urlencode({"db":"nuccore","id":acc,"rettype":"fasta","retmode":"text"})
    urls = [
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + params,
        "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?" + urllib.parse.urlencode(
            {"id":acc,"db":"nuccore","report":"fasta","retmode":"text"})
    ]
    last = None
    for url in urls:
        for attempt in range(tries):
            try:
                req = urllib.request.Request(url, headers={"User-Agent":"gc-skew-phase-coherence/primitive-closure-test contact=research"})
                with urllib.request.urlopen(req, timeout=60) as r:
                    txt = r.read().decode("utf-8")
                lines = [x.strip() for x in txt.splitlines() if x.strip()]
                if not lines or not lines[0].startswith(">"):
                    raise ValueError("response is not FASTA")
                seq = "".join(x for x in lines if not x.startswith(">")).upper()
                if seq:
                    return seq
            except Exception as e:
                last = e
                time.sleep(1.5*(attempt+1))
    raise RuntimeError(f"{acc}: failed download: {last}")

def primitive_closure(signal: np.ndarray, seed: int, n_null: int=N_NULL):
    s=np.asarray(signal,float); s=s-s.mean()
    F=np.fft.rfft(s)
    mags=np.abs(F)
    phases=np.angle(F)
    ns=np.arange(3,len(F)-2,2,dtype=int)
    w=mags[ns-2]*(mags[ns]**2)*mags[ns+2]
    sw=float(w.sum())
    if sw<=0:
        raise ValueError("zero closure weight")
    kappa=phases[ns+2]+phases[ns-2]-2.0*phases[ns]
    z1=np.exp(1j*kappa)
    z2=z1*z1
    M1=np.sum(w*z1)/sw
    M2=np.sum(w*z2)/sw
    R1=float(abs(M1))
    S1=float(np.real(M1))
    A1=float(np.angle(M1))
    R2=float(abs(M2))
    S2=float(np.real(M2))
    neff=float(sw*sw/np.sum(w*w)) if np.sum(w*w)>0 else float("nan")

    rng=np.random.default_rng(seed)
    # Exact k=1 preserved; all positive-frequency phases k>=2 below Nyquist randomized.
    # Only odd phases enter this statistic, but the full phase draw implements the original null literally.
    P=rng.uniform(0.0,2.0*np.pi,size=(n_null,len(F)))
    P[:,1]=phases[1]
    kap=P[:,ns+2]+P[:,ns-2]-2.0*P[:,ns]
    Z=np.exp(1j*kap)
    M=(Z@w)/sw
    null_R=np.abs(M)
    null_S=np.real(M)
    pR=float((1+np.sum(null_R>=R1))/(n_null+1))
    pS=float((1+np.sum(null_S>=S1))/(n_null+1))
    return {
        "M1_real":S1,
        "M1_imag":float(np.imag(M1)),
        "M1_abs":R1,
        "M1_angle_rad":A1,
        "M2_real":S2,
        "M2_imag":float(np.imag(M2)),
        "M2_abs":R2,
        "closure_quartets":int(len(ns)),
        "closure_weight_neff":neff,
        "R1_null_mean":float(np.mean(null_R)),
        "R1_null_sd":float(np.std(null_R,ddof=1)),
        "R1_null_q95":float(np.quantile(null_R,0.95)),
        "R1_null_q99":float(np.quantile(null_R,0.99)),
        "R1_p":pR,
        "S1_null_mean":float(np.mean(null_S)),
        "S1_null_sd":float(np.std(null_S,ddof=1)),
        "S1_null_q95":float(np.quantile(null_S,0.95)),
        "S1_null_q99":float(np.quantile(null_S,0.99)),
        "S1_p":pS,
    }

def bh(pvals):
    p=np.asarray(pvals,float); o=np.argsort(p); r=p[o]; m=len(r)
    q=np.minimum.accumulate((r*m/np.arange(1,m+1))[::-1])[::-1]
    out=np.empty_like(q); out[o]=np.clip(q,0,1); return out

def fisher_even(pvals):
    p=np.asarray(pvals,float)
    X=float(-2*np.log(p).sum()); df=2*len(p)
    n=df//2; y=X/2
    sf=math.exp(-y)*sum(y**j/math.factorial(j) for j in range(n))
    return X,df,sf

def frozen_J_ref(fr:dict,bins:int):
    if bins==512: return float(fr["J_max"])
    return float(fr[f"J_max_{bins}"])

def main():
    with (ROOT/"manifests"/"CONFIRMATORY_COHORT_A.csv").open(newline="") as f:
        cohort=list(csv.DictReader(f))
    with (ROOT/"results"/"confirmatory_A_results.csv").open(newline="") as f:
        frozen={r["sequence_accession"]:r for r in csv.DictReader(f)}

    rows=[]
    for i,r in enumerate(cohort):
        acc=r["sequence_accession"]; exp_len=int(r["expected_length_bp"])
        print(f"DOWNLOAD {i+1:02d}/15 {acc}",flush=True)
        seq=fetch_fasta(acc)
        if len(seq)!=exp_len:
            raise RuntimeError(f"{acc}: length {len(seq)} != expected {exp_len}")
        row={"order":int(r["order"]),"organism":r["organism"],"sequence_accession":acc,"length_bp":len(seq),
             "old_phase_p":float(frozen[acc]["k1_preserved_phase_p"])}
        for bins in BINS_LIST:
            sig=windowed_gc_skew(seq,bins=bins)
            j=float(best_projector_axis(sig)["J"])
            jref=frozen_J_ref(frozen[acc],bins)
            jerr=abs(j-jref)
            if jerr>5e-12:
                raise RuntimeError(f"{acc} bins={bins}: J mismatch {j} vs {jref}; err={jerr}")
            seed=BASE_SEED+i*100003 if bins==512 else BASE_SEED+i*100003+bins*1000003
            res=primitive_closure(sig,seed=seed)
            row[f"J_{bins}_abs_error"]=jerr
            for k,v in res.items(): row[f"{k}_{bins}"]=v
        rows.append(row)
        print("ROW",json.dumps(row,sort_keys=True),flush=True)
        time.sleep(0.35)

    summary={
        "n":len(rows),"null_reps":N_NULL,
        "primary_statistic":"R1=|sum w exp(i*kappa_n)|/sum w",
        "secondary_directional_statistic":"S1=Re[sum w exp(i*kappa_n)]/sum w",
        "kappa":"arg(c[n+2])+arg(c[n-2])-2arg(c[n]), odd n",
        "weights":"|c[n-2]|*|c[n]|^2*|c[n+2]|",
        "null":"preserve all Fourier magnitudes and exact complex k=1; randomize all higher positive-frequency phases subject to real-signal conjugacy",
        "resolution_role":"512 primary representation; 256 and 1024 robustness checks",
        "analysis_status":"post-unblinding secondary analysis of frozen cohort"
    }
    for bins in BINS_LIST:
        for stat in ("R1","S1"):
            p=[x[f"{stat}_p_{bins}"] for x in rows]
            q=bh(p)
            for x,qq in zip(rows,q): x[f"{stat}_q_BH_{bins}"]=float(qq)
            X,df,fp=fisher_even(p)
            summary[f"{bins}_{stat}"]={
                "fisher_X":X,"fisher_df":df,"fisher_p":fp,
                "raw_p_lt_0_05_n":int(sum(x<0.05 for x in p)),
                "BH_q_lt_0_05_n":int(sum(x<0.05 for x in q)),
            }
        summary[f"{bins}_descriptive"]={
            "mean_M1_abs":float(np.mean([x[f"M1_abs_{bins}"] for x in rows])),
            "median_M1_abs":float(np.median([x[f"M1_abs_{bins}"] for x in rows])),
            "mean_M1_real":float(np.mean([x[f"M1_real_{bins}"] for x in rows])),
            "median_M1_real":float(np.median([x[f"M1_real_{bins}"] for x in rows])),
            "max_J_reproduction_error":float(max(x[f"J_{bins}_abs_error"] for x in rows)),
        }

    fields=list(rows[0].keys())
    with (ROOT/"primitive_closure_results.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(rows)
    (ROOT/"primitive_closure_summary.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    print("SUMMARY",json.dumps(summary,sort_keys=True),flush=True)

if __name__=="__main__":
    main()
