#!/usr/bin/env python3
from __future__ import annotations
import csv, hashlib, json, math, sys, time, urllib.parse, urllib.request
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
from gc_symmetry_projector import windowed_gc_skew

N_NULL=4999
BASE_SEED=20261012
BINS_LIST=(256,512,1024)
IUPAC=set("ACGTRYSWKMBDHVN")
PROTOCOL_COMMIT="c3e45caad0b3a368e281ae061798b015adbc48a9"

def fetch_fasta(acc, tries=4):
    params=urllib.parse.urlencode({"db":"nuccore","id":acc,"rettype":"fasta","retmode":"text"})
    urls=[
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"+params,
        "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?"+urllib.parse.urlencode(
            {"id":acc,"db":"nuccore","report":"fasta","retmode":"text"})
    ]
    last=None
    for url in urls:
        for attempt in range(tries):
            try:
                req=urllib.request.Request(url,headers={"User-Agent":"gc-skew-phase-coherence/blinded-reserve-validation"})
                with urllib.request.urlopen(req,timeout=60) as r:
                    raw=r.read()
                txt=raw.decode("utf-8")
                lines=[x.strip() for x in txt.splitlines() if x.strip()]
                if not lines or not lines[0].startswith(">"): raise ValueError("not FASTA")
                header=lines[0][1:]
                if acc not in header and acc.split(".")[0] not in header:
                    raise ValueError(f"header mismatch: {header}")
                seq="".join(x for x in lines if not x.startswith(">")).upper()
                bad=sorted(set(seq)-IUPAC)
                if bad: raise ValueError(f"non-IUPAC symbols: {bad}")
                return raw,header,seq
            except Exception as e:
                last=e; time.sleep(1.5*(attempt+1))
    raise RuntimeError(f"{acc}: download failed: {last}")

def primitive_stats(signal, seed, n_null=N_NULL):
    s=np.asarray(signal,float); s=s-s.mean()
    F=np.fft.rfft(s)
    mags=np.abs(F); phases=np.angle(F)
    ns=np.arange(3,len(F)-2,2,dtype=int)
    w=mags[ns-2]*(mags[ns]**2)*mags[ns+2]
    sw=float(w.sum())
    if sw<=0: raise ValueError("zero closure weight")
    kappa=phases[ns+2]+phases[ns-2]-2.0*phases[ns]
    z1=np.exp(1j*kappa); z2=z1*z1
    M1=np.sum(w*z1)/sw; M2=np.sum(w*z2)/sw
    S1=float(np.real(M1)); R1=float(abs(M1))
    neff=float(sw*sw/np.sum(w*w)) if np.sum(w*w)>0 else float("nan")

    rng=np.random.default_rng(seed)
    P=rng.uniform(0.0,2.0*np.pi,size=(n_null,len(F)))
    P[:,1]=phases[1]
    kap=P[:,ns+2]+P[:,ns-2]-2.0*P[:,ns]
    M=(np.exp(1j*kap)@w)/sw
    nullS=np.real(M)
    p=float((1+np.sum(nullS>=S1))/(n_null+1))
    return {
        "S1":S1,"R1":R1,"M1_imag":float(np.imag(M1)),"M1_angle_rad":float(np.angle(M1)),
        "M2_real":float(np.real(M2)),"M2_abs":float(abs(M2)),
        "closure_quartets":int(len(ns)),"closure_weight_neff":neff,
        "null_mean_S1":float(np.mean(nullS)),"null_sd_S1":float(np.std(nullS,ddof=1)),
        "null_q95_S1":float(np.quantile(nullS,0.95)),"null_q99_S1":float(np.quantile(nullS,0.99)),
        "S1_p":p
    }

def bh(p):
    p=np.asarray(p,float); o=np.argsort(p); r=p[o]; m=len(r)
    q=np.minimum.accumulate((r*m/np.arange(1,m+1))[::-1])[::-1]
    out=np.empty_like(q); out[o]=np.clip(q,0,1); return out

def chi2_sf_even(x,df):
    n=df//2; y=x/2
    return math.exp(-y)*sum(y**j/math.factorial(j) for j in range(n))

def main():
    with (ROOT/"manifests"/"CONFIRMATORY_RESERVE_QUEUE.csv").open(newline="") as f:
        cohort=list(csv.DictReader(f))
    if len(cohort)!=5: raise RuntimeError(f"reserve cohort size changed: {len(cohort)}")

    rows=[]
    for r in cohort:
        order=int(r["order"]); acc=r["sequence_accession"]; expected=int(r["expected_length_bp"])
        print(f"OPENING BLINDED RESERVE {order}/5 {acc}",flush=True)
        raw,header,seq=fetch_fasta(acc)
        if len(seq)!=expected: raise RuntimeError(f"{acc}: length {len(seq)} != {expected}")
        row={
            "order":order,"organism":r["organism"],"sequence_accession":acc,
            "expected_length_bp":expected,"fasta_header":header,
            "raw_sha256":hashlib.sha256(raw).hexdigest(),
            "canonical_sha256":hashlib.sha256(seq.encode()).hexdigest(),
            "sequence_length":len(seq),"provenance_status":"PASS"
        }
        for bins in BINS_LIST:
            sig=windowed_gc_skew(seq,bins=bins)
            seed=BASE_SEED+100003*order
            if bins==256: seed+=256000003
            if bins==1024: seed+=1024000003
            res=primitive_stats(sig,seed)
            for k,v in res.items(): row[f"{k}_{bins}"]=v
        rows.append(row)
        print("RESULT",json.dumps(row,sort_keys=True),flush=True)
        time.sleep(0.35)

    p512=[x["S1_p_512"] for x in rows]
    q512=bh(p512)
    for x,q in zip(rows,q512): x["S1_q_BH_512"]=float(q)
    X=float(-2*np.log(np.asarray(p512)).sum()); fp=chi2_sf_even(X,10)
    summary={
        "protocol_commit":PROTOCOL_COMMIT,"cohort_n":5,"primary_resolution":512,
        "primary_statistic":"S1=Re(sum w exp(i*kappa_n)/sum w)",
        "null_reps":N_NULL,"per_genome_p":[float(x) for x in p512],
        "fisher_X":X,"fisher_df":10,"fisher_p":fp,
        "primary_success_threshold":0.05,"BLINDED_VALIDATION_POSITIVE":bool(fp<0.05),
        "raw_p_lt_0_05_n":int(sum(x<0.05 for x in p512)),
        "BH_q_lt_0_05_n":int(sum(x<0.05 for x in q512)),
        "mean_S1_512":float(np.mean([x["S1_512"] for x in rows])),
        "median_S1_512":float(np.median([x["S1_512"] for x in rows])),
        "all_provenance_pass":all(x["provenance_status"]=="PASS" for x in rows),
        "robustness":{}
    }
    for bins in (256,1024):
        pp=[x[f"S1_p_{bins}"] for x in rows]
        XX=float(-2*np.log(np.asarray(pp)).sum())
        summary["robustness"][str(bins)]={
            "fisher_X":XX,"fisher_df":10,"fisher_p":chi2_sf_even(XX,10),
            "raw_p_lt_0_05_n":int(sum(x<0.05 for x in pp)),
            "mean_S1":float(np.mean([x[f"S1_{bins}"] for x in rows]))
        }

    with (ROOT/"primitive_closure_results.csv").open("w",newline="") as f:
        fields=list(rows[0].keys()); w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(rows)
    (ROOT/"primitive_closure_summary.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    print("PRIMARY_SUMMARY",json.dumps(summary,sort_keys=True),flush=True)

if __name__=="__main__":
    main()
