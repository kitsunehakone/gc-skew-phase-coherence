#!/usr/bin/env python3
from __future__ import annotations
import argparse, hashlib, math, re
from pathlib import Path
import numpy as np
import pandas as pd
from gc_symmetry_projector import windowed_gc_skew
from skewdb_projector_benchmark import analyze_signal

IUPAC=set("ACGTRYSWKMBDHVN")

def read_fasta(path):
    raw=Path(path).read_bytes()
    text=raw.decode("utf-8")
    heads=[]; parts=[]
    for line in text.splitlines():
        s=line.strip()
        if not s: continue
        if s.startswith(">"): heads.append(s[1:])
        else: parts.append(re.sub(r"\s+","",s).upper())
    if not heads or not parts: raise ValueError("not FASTA")
    seq="".join(parts)
    bad=sorted(set(seq)-IUPAC)
    if bad: raise ValueError(f"non-IUPAC symbols: {bad}")
    return raw,heads,seq

def sha(x): return hashlib.sha256(x).hexdigest()

def locate(d,acc):
    d=Path(d)
    for p in list(d.glob(f"{acc}*"))+list(d.glob("*")):
        if not p.is_file(): continue
        try:
            _,h,_=read_fasta(p)
            if any(acc in x or acc.split(".")[0] in x.split()[0] for x in h):
                return p
        except Exception: pass
    raise FileNotFoundError(acc)

def bh(p):
    p=np.asarray(p,float); o=np.argsort(p); r=p[o]; m=len(r)
    q=np.minimum.accumulate((r*m/np.arange(1,m+1))[::-1])[::-1]
    out=np.empty_like(q); out[o]=np.clip(q,0,1); return out

def wilson(k,n,z=1.959963984540054):
    ph=k/n; d=1+z*z/n
    c=(ph+z*z/(2*n))/d
    h=z*math.sqrt(ph*(1-ph)/n+z*z/(4*n*n))/d
    return max(0,c-h),min(1,c+h)

def chi2_sf_even(x,df):
    n=df//2; y=x/2
    return math.exp(-y)*sum(y**j/math.factorial(j) for j in range(n))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--cohort",default="CONFIRMATORY_COHORT_A.csv")
    ap.add_argument("--fasta-dir",required=True)
    ap.add_argument("--outdir",required=True)
    ap.add_argument("--null-reps",type=int,default=4999)
    ap.add_argument("--seed",type=int,default=20260920)
    a=ap.parse_args()

    cohort=pd.read_csv(a.cohort); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    rows=[]; prov=[]
    for i,r in cohort.iterrows():
        acc=str(r.sequence_accession); p=locate(a.fasta_dir,acc)
        raw,heads,seq=read_fasta(p)
        if len(seq)!=int(r.expected_length_bp):
            raise ValueError(f"{acc}: length {len(seq)} != {int(r.expected_length_bp)}")
        prov.append(dict(sequence_accession=acc,file=p.name,length_bp=len(seq),
                         raw_sha256=sha(raw),canonical_sha256=sha(seq.encode()),status="PASS"))
        row=dict(order=int(r["order"]),organism=r.organism,sequence_accession=acc)
        sig=windowed_gc_skew(seq,bins=512)
        row.update(analyze_signal(sig,n_null=a.null_reps,seed=a.seed+i*100003))
        a512=row["projector_axis_fraction"]
        for M in (256,1024):
            rr=analyze_signal(windowed_gc_skew(seq,bins=M),n_null=0)
            d=abs(rr["projector_axis_fraction"]-a512)%0.5; d=min(d,0.5-d)
            row[f"axis_shift_{M}_degrees"]=d*360
            row[f"J_max_{M}"]=rr["J_max"]
        row["resolution_stable"]=(row["axis_shift_256_degrees"]<=2*(360/512)
                                  and row["axis_shift_1024_degrees"]<=2*(360/1024))
        rows.append(row)

    df=pd.DataFrame(rows)
    df["phase_q_BH"]=bh(df["k1_preserved_phase_p"])
    df.to_csv(out/"confirmatory_A_results.csv",index=False)
    pd.DataFrame(prov).to_csv(out/"confirmatory_A_provenance.csv",index=False)

    p=df["k1_preserved_phase_p"].to_numpy(float)
    X=float(-2*np.log(p).sum()); fp=chi2_sf_even(X,2*len(p))
    k=int((p<0.05).sum()); lo,hi=wilson(k,len(p))
    summary=pd.DataFrame([dict(
        n=len(p),fisher_statistic=X,fisher_df=2*len(p),fisher_p=fp,
        confirmatory_positive=(fp<0.05),raw_p_lt_0_05_n=k,
        raw_p_lt_0_05_fraction=k/len(p),wilson95_lo=lo,wilson95_hi=hi,
        BH_q_lt_0_05_n=int((df.phase_q_BH<0.05).sum()),
        resolution_stable_n=int(df.resolution_stable.sum())
    )])
    summary.to_csv(out/"confirmatory_A_summary.csv",index=False)
    print(summary.to_string(index=False))

if __name__=="__main__": main()
