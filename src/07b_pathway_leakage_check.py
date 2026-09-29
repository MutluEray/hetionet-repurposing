"""Day 6: do disease-named pathways make the pathway metapath circular?

Path graphs showed pathways such as 'Alzheimers Disease' and 'Parkinsons Disease Pathway' among the
evidence for Alzheimer's/Parkinson's candidates. Genes in a pathway named after the disease are
close to the disease genes by construction, so paths through them add little independent evidence.
This script drops every pathway whose name contains a Hetionet disease name and compares:

  * the pan-disease benchmark (CbGpPWpGaD and the binding aggregate, all pathways vs no disease pathways)
  * the AD / PD top-15 candidate lists (DWPC, binding aggregate)

Usage:
    python3 src/07b_pathway_leakage_check.py
    python3 src/07b_pathway_leakage_check.py --w 0.6
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from dwpc import aggregate_vector, build_matrices, compute_dwpc, coverage, rank_disease
from evalutils import metrics, rho, summarize
from utils import OUT_TAB, find_nodes, load_processed

BIND = ["CbGaD", "CbGiGaD", "CbGpPWpGaD"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--w", type=float, default=0.6)
    ap.add_argument("--min-pos", type=int, default=3)
    ap.add_argument("--min-nonzero", type=int, default=20)
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    nodes, edges = load_processed()
    variants = {"all": build_matrices(nodes, edges),
                "nodis": build_matrices(nodes, edges, drop_disease_pathways=True)}
    dropped = variants["nodis"].dropped_pathways
    n_all = variants["all"].A["GpPW"].nnz
    n_left = variants["nodis"].A["GpPW"].nnz
    print(f"Disease-named pathways removed: {len(dropped)} of {len(variants['all'].tables['Pathway'])} "
          f"({n_all - n_left:,} of {n_all:,} gene-pathway edges)")
    print("Examples:", "; ".join(dropped[:40]))
    key = [p for p in dropped if "alzheimer" in p.lower() or "parkinson" in p.lower()]
    print("Alzheimer/Parkinson-named pathways:", key)

    A0 = variants["all"].A
    deg = np.asarray(A0["CbG"].sum(axis=1)).ravel().astype(float)
    CtD = A0["CtD"].toarray() > 0
    bench = np.where(CtD.sum(axis=0) >= args.min_pos)[0]
    dis = variants["all"].tables["Disease"]

    scores, covs = {}, {}
    for name, m in variants.items():
        scores[name] = compute_dwpc(m.A, args.w)
        covs[name] = coverage(m.A)

    rows = []
    for d in bench:
        pos = CtD[:, d]
        vecs = {"degree_CbG": deg}
        for name in variants:
            sc, cv = scores[name], covs[name]
            v = sc["CbGpPWpGaD"][:, d].astype(float).copy()
            v[~cv["CbGpPWpGaD"]] = np.nan
            vecs[f"CbGpPWpGaD|{name}"] = v
            vecs[f"agg_bind|{name}"] = aggregate_vector(sc, cv, d, BIND, "mean", None, args.min_nonzero)[0]
        for mname, v in vecs.items():
            auc, med, h5, h10 = metrics(v, pos)
            rows.append(dict(disease_id=dis.loc[d, "id"], disease=dis.loc[d, "name"], n_pos=int(pos.sum()),
                             method=mname, auroc=auc, median_rank_pct=med, hit5=h5, hit10=h10,
                             rho_degree=rho(v, deg)))
    per = pd.DataFrame(rows)
    summ = summarize(per, args.boot, rng)
    summ.to_csv(OUT_TAB / "eval_pathway_leakage.csv", index=False)
    print(f"\n=== Pan-disease benchmark, w={args.w:g}, {len(bench)} diseases ===")
    print(summ[["method", "mean_auroc", "ci_lo", "ci_hi", "wins_vs_degree", "p_vs_degree",
                "rho_degree"]].round(3).to_string(index=False))

    piv = per.pivot(index="disease_id", columns="method", values="auroc").fillna(0.5)
    for m in ("CbGpPWpGaD", "agg_bind"):
        diff = piv[f"{m}|nodis"] - piv[f"{m}|all"]
        try:
            p = wilcoxon(diff)[1]
        except Exception:
            p = float("nan")
        print(f"Paired change in AUROC for {m} when disease pathways are dropped: "
              f"mean {diff.mean():+.4f}, wins {int((diff > 0).sum())}/{len(diff)}, Wilcoxon p={p:.3f}")

    # ---------------- AD / PD candidate lists ----------------
    for q in ("alzheimer", "parkinson"):
        hits = find_nodes(nodes, q, kind="Disease")
        if hits.empty:
            continue
        d = variants["all"].local_index(int(hits.iloc[0]["idx"]))
        tops = {}
        for name, m in variants.items():
            df, _ = rank_disease(m, scores[name], covs[name], d, BIND, mode="mean", min_nonzero=args.min_nonzero)
            tops[name] = df[~df["known"] & df["score"].notna()].head(15)["compound"].tolist()
        print(f"\n=== {hits.iloc[0]['name']}: top-15 (DWPC binding aggregate, known drugs excluded) ===")
        print(pd.DataFrame({"all pathways": tops["all"], "no disease pathways": tops["nodis"]}).to_string())
        print(f"Overlap: {len(set(tops['all']) & set(tops['nodis']))}/15")


if __name__ == "__main__":
    main()
