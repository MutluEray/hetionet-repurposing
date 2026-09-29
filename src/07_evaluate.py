"""Day 4: evaluation, baselines, damping sensitivity and hub audit.

Question this script answers: does the graph-based score add anything beyond
"this compound binds many genes" (promiscuity)?

Pan-disease benchmark
    Positives = Hetionet 'treats' (CtD) edges; diseases with >= --min-pos positives.
    Negatives = all other compounds (UNLABELED, so some are true positives: AUROC is
    conservative). Compounds with no data for a method are ranked at the bottom.
    Metrics per disease: AUROC, median rank percentile of positives, hit@5% / hit@10%.
    Macro-averaged over diseases with a bootstrap CI (resampling diseases), a paired
    Wilcoxon test and win-rate versus the degree baseline.
    Leakage note: a positive's own CtD edge is never used to score it (CrC / DrD have no
    self-loops), but degree weights include it; CrCtD/CtDrD use OTHER known treatments.
    Circularity note: drug targets (CbG) and disease genes (DaG) come from overlapping
    literature, so target-based scores partly encode the indications.

Baselines: degree (# genes the compound binds), raw (w=0) metapath scores.
Also: damping sweep (macro AUROC and Spearman with degree vs w), AD/PD evaluation with
positives = treats U palliates (few positives -> descriptive only), and a hub audit.

Outputs (outputs/tables/): eval_pan_disease_summary.csv, eval_pan_disease_per_disease.csv,
eval_damping.csv, eval_ad_pd.csv

Usage:
    python3 src/07_evaluate.py
    python3 src/07_evaluate.py --min-pos 3 --w 0.4 --boot 2000
"""
from __future__ import annotations

import argparse
import time

import numpy as np
import pandas as pd

from dwpc import (METAPATHS, MECHANISTIC, aggregate_vector, build_matrices,
                  compute_dwpc, coverage, known_mask)
from evalutils import metrics, rho, summarize
from utils import OUT_TAB, find_nodes, load_processed

AGG_VARIANTS = {
    "agg_mech_mean": (MECHANISTIC, "mean"),
    "agg_mech_groups": (MECHANISTIC, "groups"),
    "agg_all_mean": (METAPATHS, "mean"),
    "agg_all_groups": (METAPATHS, "groups"),
}


def vectors(scores, cov, d, suffix, min_nonzero):
    out = {}
    for mp in METAPATHS:
        v = scores[mp][:, d].astype(float).copy()
        v[~cov[mp]] = np.nan
        out[f"{mp}@{suffix}"] = v
    for name, (mps, mode) in AGG_VARIANTS.items():
        s, _, _ = aggregate_vector(scores, cov, d, mps, mode, None, min_nonzero)
        out[f"{name}@{suffix}"] = s
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--w", type=float, default=0.4)
    ap.add_argument("--min-pos", type=int, default=3)
    ap.add_argument("--min-nonzero", type=int, default=20)
    ap.add_argument("--w-grid", default="0,0.2,0.4,0.6,0.8,1.0")
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    t0 = time.time()
    rng = np.random.default_rng(args.seed)

    nodes, edges = load_processed()
    mats = build_matrices(nodes, edges)
    A = mats.A
    cov = coverage(A)
    dis = mats.tables["Disease"]
    deg = np.asarray(A["CbG"].sum(axis=1)).ravel().astype(float)
    CtD = A["CtD"].toarray() > 0
    bench = np.where(CtD.sum(axis=0) >= args.min_pos)[0]
    print(f"Pan-disease benchmark: {len(bench)} diseases with >= {args.min_pos} treats edges "
          f"({int(CtD[:, bench].sum())} positives)")

    main_tag = f"w{args.w:g}"
    S_main = compute_dwpc(A, args.w)
    S_raw = compute_dwpc(A, 0.0)

    # ---------------- pan-disease benchmark ----------------
    rows = []
    for d in bench:
        pos = CtD[:, d]
        vecs = {"degree_CbG": deg}
        vecs.update(vectors(S_raw, cov, d, "raw", args.min_nonzero))
        vecs.update(vectors(S_main, cov, d, main_tag, args.min_nonzero))
        for name, v in vecs.items():
            auc, med, h5, h10 = metrics(v, pos)
            rows.append(dict(disease_id=dis.loc[d, "id"], disease=dis.loc[d, "name"], n_pos=int(pos.sum()),
                             method=name, auroc=auc, median_rank_pct=med, hit5=h5, hit10=h10,
                             rho_degree=rho(v, deg)))
    per = pd.DataFrame(rows)
    per.to_csv(OUT_TAB / "eval_pan_disease_per_disease.csv", index=False)

    summ = summarize(per, args.boot, rng)
    summ.to_csv(OUT_TAB / "eval_pan_disease_summary.csv", index=False)
    print("\n=== Pan-disease benchmark (macro over diseases; 0.5 = chance) ===")
    print(summ.round(3).to_string(index=False))

    # ---------------- damping sweep ----------------
    grid = [float(x) for x in args.w_grid.split(",")]
    keep = ["CbGaD", "CbGiGaD", "CbGpPWpGaD", "agg_mech_mean", "agg_mech_groups"]
    drows = []
    for wg in grid:
        sc = compute_dwpc(A, wg)
        acc = {k: [] for k in keep}
        rr = {k: [] for k in keep}
        for d in bench:
            pos = CtD[:, d]
            vecs = vectors(sc, cov, d, "x", args.min_nonzero)
            for k in keep:
                v = vecs[f"{k}@x"]
                acc[k].append(np.nan_to_num(metrics(v, pos)[0], nan=0.5))
                rr[k].append(rho(v, deg))
        for k in keep:
            drows.append(dict(w=wg, method=k, mean_auroc=np.mean(acc[k]), mean_rho_degree=np.nanmean(rr[k])))
    dm = pd.DataFrame(drows)
    dm.to_csv(OUT_TAB / "eval_damping.csv", index=False)
    print("\n=== Damping sweep: mean AUROC ===")
    print(dm.pivot(index="w", columns="method", values="mean_auroc").round(3).to_string())
    print("\n=== Damping sweep: mean Spearman(score, # bound genes) ===")
    print(dm.pivot(index="w", columns="method", values="mean_rho_degree").round(3).to_string())

    # ---------------- AD / PD ----------------
    ad_rows, hub_rows = [], []
    for q in ("alzheimer", "parkinson"):
        hits = find_nodes(nodes, q, kind="Disease")
        if hits.empty:
            continue
        d = mats.local_index(int(hits.iloc[0]["idx"]))
        pos = known_mask(A, d)
        vecs = {"degree_CbG": deg}
        vecs.update(vectors(S_raw, cov, d, "raw", args.min_nonzero))
        vecs.update(vectors(S_main, cov, d, main_tag, args.min_nonzero))
        for name, v in vecs.items():
            auc, med, h5, h10 = metrics(v, pos)
            ad_rows.append(dict(disease=hits.iloc[0]["name"], n_pos=int(pos.sum()), method=name,
                                auroc=auc, median_rank_pct=med, hit5=h5, hit10=h10))
        for name in ("agg_mech_mean", "agg_mech_groups"):
            v = vecs[f"{name}@{main_tag}"]
            order = np.argsort(-np.nan_to_num(v, nan=-1))
            top = [i for i in order if not pos[i] and np.isfinite(v[i])][:20]
            hub_rows.append(dict(disease=hits.iloc[0]["name"], method=name,
                                 median_degree_top20=float(np.median(deg[top])),
                                 median_degree_all_covered=float(np.median(deg[cov["CbGaD"]]))))
    adpd = pd.DataFrame(ad_rows)
    adpd.to_csv(OUT_TAB / "eval_ad_pd.csv", index=False)
    print("\n=== AD / PD (positives = treats U palliates; few positives -> descriptive) ===")
    for dname, g in adpd.groupby("disease"):
        print(f"\n{dname}  (n_pos={g['n_pos'].iloc[0]})")
        print(g.drop(columns=["disease", "n_pos"]).sort_values("auroc", ascending=False)
              .round(3).to_string(index=False))
    print("\n=== Hub audit (top-20 candidates, known drugs excluded) ===")
    print(pd.DataFrame(hub_rows).round(1).to_string(index=False))
    print(f"\nDone in {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
