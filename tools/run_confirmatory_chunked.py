#!/usr/bin/env python3
"""Additive execution wrapper for the frozen confirmatory runner.

IMPORTANT PROVENANCE STATUS
---------------------------
This convenience wrapper was created *after* Confirmatory Cohort A was opened.
It does not change the frozen statistic, null model, cohort, seed schedule, or
scientific code.  It only invokes the frozen pre-unblinding runner one genome
at a time so jobs can fit environments with short wall-clock limits, then
recombines the resulting per-genome outputs.

The canonical scientific implementation remains:
    code/run_confirmatory_prevalence.py
    code/gc_symmetry_projector.py
    code/skewdb_projector_benchmark.py

For original cohort row i (zero-based), the frozen runner used seed
    base_seed + i * 100003.
When a single-row cohort is invoked, this wrapper supplies that adjusted seed
so every chromosome receives the same RNG stream as in the original full run.
"""
from __future__ import annotations
import argparse, math, subprocess, sys, tempfile
from pathlib import Path
import numpy as np
import pandas as pd


def bh_adjust(p):
    p=np.asarray(p,float); o=np.argsort(p); r=p[o]; m=len(r)
    q=np.minimum.accumulate((r*m/np.arange(1,m+1))[::-1])[::-1]
    out=np.empty_like(q); out[o]=np.clip(q,0,1); return out


def wilson(k,n,z=1.959963984540054):
    ph=k/n; d=1+z*z/n
    c=(ph+z*z/(2*n))/d
    h=z*math.sqrt(ph*(1-ph)/n+z*z/(4*n*n))/d
    return max(0,c-h),min(1,c+h)


def chi2_sf_even(x,df):
    if not isinstance(df, (int, np.integer)) or df <= 0 or df % 2:
        raise ValueError(f"chi2_sf_even requires a positive even integer df; got {df!r}")
    n=df//2; y=x/2
    return math.exp(-y)*sum(y**j/math.factorial(j) for j in range(n))


def frozen_runner(repo: Path) -> Path:
    p=repo/'code'/'run_confirmatory_prevalence.py'
    if not p.exists(): raise FileNotFoundError(p)
    return p


def run_one(args):
    repo=Path(args.repo).resolve()
    cohort=pd.read_csv(args.cohort)
    matches=cohort.index[cohort['sequence_accession'].astype(str)==args.genome].tolist()
    if len(matches)!=1:
        raise SystemExit(f'Expected exactly one cohort row for {args.genome}; found {len(matches)}')
    idx=matches[0]
    row=cohort.iloc[[idx]]
    adjusted_seed=args.seed + idx*100003
    out=Path(args.outdir).resolve()/'per_genome'/args.genome
    out.mkdir(parents=True,exist_ok=True)
    single=out/'single_row_cohort.csv'
    row.to_csv(single,index=False)
    cmd=[sys.executable,str(frozen_runner(repo)),
         '--cohort',str(single),'--fasta-dir',str(Path(args.fasta_dir).resolve()),
         '--outdir',str(out),'--null-reps',str(args.null_reps),'--seed',str(adjusted_seed)]
    print('Running:', ' '.join(cmd), flush=True)
    subprocess.run(cmd,check=True)
    meta=pd.DataFrame([{
        'sequence_accession':args.genome,
        'original_zero_based_index':idx,
        'base_seed':args.seed,
        'adjusted_seed':adjusted_seed,
        'null_reps':args.null_reps,
    }])
    meta.to_csv(out/'wrapper_execution_metadata.csv',index=False)
    print(f'Completed {args.genome} with adjusted seed {adjusted_seed}')


def merge(args):
    cohort=pd.read_csv(args.cohort).reset_index(drop=True)
    root=Path(args.outdir).resolve()/'per_genome'
    rows=[]; prov=[]
    for idx,crow in cohort.iterrows():
        acc=str(crow.sequence_accession)
        d=root/acc
        rp=d/'confirmatory_A_results.csv'; pp=d/'confirmatory_A_provenance.csv'
        if not rp.exists():
            raise SystemExit(f'Missing result for {acc}: {rp}')
        if not pp.exists():
            raise SystemExit(f'Missing provenance file for {acc}: {pp}')
        r=pd.read_csv(rp)
        if len(r)!=1:
            raise SystemExit(f'Expected exactly one result row for {acc}; found {len(r)}')
        required = {
            'sequence_accession', 'J_max', 'k1_preserved_phase_p',
            'projector_axis_fraction', 'first_harmonic_power_fraction',
            'higher_odd_energy_fraction', 'odd_axis_coherence',
            'resolution_stable'
        }
        missing = sorted(required - set(r.columns))
        if missing:
            raise SystemExit(f'Truncated/incomplete result for {acc}; missing columns: {missing}')
        stored_acc = str(r.iloc[0]['sequence_accession'])
        if stored_acc != acc:
            raise SystemExit(f'Result accession mismatch: expected {acc}, found {stored_acc}')
        numeric_required = [
            'J_max', 'k1_preserved_phase_p', 'projector_axis_fraction',
            'first_harmonic_power_fraction', 'higher_odd_energy_fraction',
            'odd_axis_coherence'
        ]
        for col in numeric_required:
            val = pd.to_numeric(r.iloc[0][col], errors='coerce')
            if not np.isfinite(val):
                raise SystemExit(f'Truncated/incomplete result for {acc}; {col} is not finite')
        pval = float(r.iloc[0]['k1_preserved_phase_p'])
        pmin = 1.0/(args.null_reps+1.0)
        if not (pmin <= pval <= 1.0):
            raise SystemExit(
                f'Invalid Monte Carlo p-value for {acc}: {pval}; expected [{pmin}, 1]'
            )
        pr=pd.read_csv(pp)
        if len(pr)!=1:
            raise SystemExit(f'Expected exactly one provenance row for {acc}; found {len(pr)}')
        if 'sequence_accession' not in pr.columns or str(pr.iloc[0]['sequence_accession']) != acc:
            raise SystemExit(f'Provenance accession mismatch for {acc}')
        rec=r.iloc[0].to_dict()
        rec.update({
            'order':int(crow['order']),
            'organism':crow['organism'],
            'genus':crow.get('genus',''),
            'sequence_accession':acc,
            'expected_length_bp':int(crow['expected_length_bp']),
            'analysis_bins':512,
        })
        # Remove the one-row BH adjustment emitted by the frozen runner;
        # cohort-wide BH is recomputed below.
        rec.pop('phase_q_BH',None)
        rows.append(rec)
        x=pr.iloc[0].to_dict(); x['order']=int(crow['order']); prov.append(x)

    df=pd.DataFrame(rows).sort_values('order').reset_index(drop=True)
    df['k1_preserved_phase_q_BH']=bh_adjust(df['k1_preserved_phase_p'].to_numpy(float))
    if 'k1_preserved_null_q95_J' in df.columns:
        df['J_minus_null_q95']=df['J_max']-df['k1_preserved_null_q95_J']
    p=df['k1_preserved_phase_p'].to_numpy(float)
    X=float(-2*np.log(p).sum()); fp=chi2_sf_even(X,2*len(p))
    k=int((p<0.05).sum()); lo,hi=wilson(k,len(p))
    summary=pd.DataFrame([dict(
        n=len(p),null_replicates_per_genome=args.null_reps,
        fisher_statistic=X,fisher_df=2*len(p),fisher_p=fp,
        confirmatory_positive=(fp<0.05),raw_p_lt_0_05_n=k,
        raw_p_lt_0_05_fraction=k/len(p),wilson95_lo=lo,wilson95_hi=hi,
        BH_q_lt_0_05_n=int((df.k1_preserved_phase_q_BH<0.05).sum()),
        resolution_stable_n=int(df.resolution_stable.sum()),
        median_J_max=float(df.J_max.median()),
        median_P1=float(df.first_harmonic_power_fraction.median()),
        median_higher_odd=float(df.higher_odd_energy_fraction.median()),
        median_C_odd=float(df.odd_axis_coherence.median()),
    )])
    final=Path(args.outdir).resolve(); final.mkdir(parents=True,exist_ok=True)
    df.to_csv(final/'confirmatory_A_results_chunked.csv',index=False)
    summary.to_csv(final/'confirmatory_A_summary_chunked.csv',index=False)
    if prov: pd.DataFrame(prov).sort_values('order').to_csv(final/'confirmatory_A_provenance_chunked.csv',index=False)
    print(summary.to_string(index=False))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--repo',default=Path(__file__).resolve().parents[1])
    ap.add_argument('--cohort',required=True)
    ap.add_argument('--fasta-dir')
    ap.add_argument('--outdir',required=True)
    ap.add_argument('--null-reps',type=int,default=4999)
    ap.add_argument('--seed',type=int,default=20260920)
    g=ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--genome',help='Run one frozen cohort accession using its original seed offset.')
    g.add_argument('--merge',action='store_true',help='Merge all completed per-genome jobs and recompute cohort statistics.')
    a=ap.parse_args()
    if a.genome:
        if not a.fasta_dir: ap.error('--fasta-dir is required with --genome')
        run_one(a)
    else:
        merge(a)

if __name__=='__main__': main()
