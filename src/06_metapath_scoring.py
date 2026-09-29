"""Day 3 (v2): metapath DWPC scoring and explainable ranking for AD / PD.

Builds sparse matrices from the full Hetionet edge list (seconds, tiny memory),
computes degree-weighted path counts per metapath for ALL compound x disease pairs,
aggregates them into a coverage-aware percentile-rank score, and writes:

    data/processed/dwpc_w<w>.npz                 all metapath matrices
    outputs/tables/ranking_<slug>_<tag>.csv      candidates (known drugs excluded) + per-metapath columns
    outputs/tables/known_ranks_<slug>_<tag>.csv  where the graph's known drugs land in the ranking

It also prints raw path counts (w=0) next to DWPC (w=<w>) so the hub effect is visible.
Metapaths with fewer than --min-nonzero non-zero compounds for a disease are dropped for it.

Usage:
    python3 src/06_metapath_scoring.py
    python3 src/06_metapath_scoring.py --mode groups
    python3 src/06_metapath_scoring.py --exclude CrCtD,CtDrD
"""
from __future__ import annotations

import argparse
import time

import numpy as np
import pandas as pd

from dwpc import METAPATHS, build_matrices, compute_dwpc, coverage, rank_disease
from utils import DATA_PROC, OUT_TAB, find_nodes, load_processed


def show(df: pd.DataFrame, used: list, n: int) -> None:
    cols = ["rank", "compound", "score", "n_covered"] + [f"pct_{m}" for m in used]
    out = df[~df["known"] & df["score"].notna()].head(n)[cols].copy()
    out["rank"] = out["rank"].astype(int)
    print(out.round(3).to_string(index=False))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--w", type=float, default=0.4, help="DWPC damping exponent (Rephetio used ~0.4)")
    ap.add_argument("--diseases", default="alzheimer,parkinson")
    ap.add_argument("--exclude", default="", help="comma-separated metapaths to leave out")
    ap.add_argument("--mode", default="mean", choices=["mean", "groups"])
    ap.add_argument("--min-nonzero", type=int, default=20)
    ap.add_argument("--top", type=int, default=25)
    args = ap.parse_args()

    excluded = {x.strip() for x in args.exclude.split(",") if x.strip()}
    mps = [m for m in METAPATHS if m not in excluded]
    tag = f"w{args.w:g}_{args.mode}" + ("_excl-" + "-".join(sorted(excluded)) if excluded else "")

    t0 = time.time()
    nodes, edges = load_processed()
    mats = build_matrices(nodes, edges)
    scores = compute_dwpc(mats.A, args.w)
    raw = compute_dwpc(mats.A, 0.0)
    cov = coverage(mats.A)
    print(f"Matrices + DWPC (w={args.w:g} and raw) computed in {time.time() - t0:.1f}s")
    print(f"Metapaths considered: {', '.join(mps)}   mode={args.mode}")

    comp = mats.tables["Compound"]
    dis = mats.tables["Disease"]
    np.savez_compressed(DATA_PROC / f"dwpc_w{args.w:g}.npz", **scores,
                        compound_ids=comp["id"].to_numpy(), disease_ids=dis["id"].to_numpy())
    deg = np.asarray(mats.A["CbG"].sum(axis=1)).ravel()

    for query in [q.strip() for q in args.diseases.split(",") if q.strip()]:
        hits = find_nodes(nodes, query, kind="Disease")
        if hits.empty:
            print(f"\nNo disease matches '{query}'")
            continue
        row = hits.iloc[0]
        d = mats.local_index(int(row["idx"]))
        print(f"\n{'=' * 78}\n{row['name']}  ({row['id']})\n{'=' * 78}")

        full, used = rank_disease(mats, scores, cov, d, mps, mode=args.mode, min_nonzero=args.min_nonzero)
        full_raw, used_raw = rank_disease(mats, raw, cov, d, mps, mode=args.mode, min_nonzero=args.min_nonzero)
        dropped = [m for m in mps if m not in used]
        print(f"Metapaths used: {used}   dropped (too few non-zero compounds): {dropped}")

        print(f"\nTop {args.top} candidates, DWPC w={args.w:g} (known drugs excluded):")
        show(full, used, args.top)

        print("\nTop 15 candidates with RAW path counts (w=0) for contrast:")
        show(full_raw, used_raw, 15)
        top = full[~full["known"] & full["score"].notna()].head(20)
        top_raw = full_raw[~full_raw["known"] & full_raw["score"].notna()].head(20)
        print(f"\nTop-20 overlap raw vs DWPC: {len(set(top['compound']) & set(top_raw['compound']))}/20")
        d_all = deg[cov["CbGaD"]]
        d_top = deg[top["cidx"].to_numpy()]
        print(f"Hub audit: median # bound genes - top-20 candidates {np.median(d_top):.0f} "
              f"vs all covered compounds {np.median(d_all):.0f}")

        slug = query.lower().replace(" ", "_")
        full[~full["known"] & full["score"].notna()].head(200).to_csv(
            OUT_TAB / f"ranking_{slug}_{tag}.csv", index=False)

        ranked = full[full["score"].notna()]
        known = full[full["known"] & full["score"].notna()].copy()
        known["n_ranked"] = len(ranked)
        known["rank_pct"] = known["rank"] / len(ranked)
        keep = ["compound", "rank", "n_ranked", "rank_pct", "score", "n_covered"] + [f"pct_{m}" for m in used]
        known[keep].to_csv(OUT_TAB / f"known_ranks_{slug}_{tag}.csv", index=False)
        print(f"\nKnown drugs (treats/palliates) in this ranking - {len(ranked)} compounds ranked:")
        print(known[keep].sort_values("rank").round(3).to_string(index=False))
        if len(known):
            print(f"Median rank percentile of known drugs (lower = better; random ~0.5): "
                  f"{known['rank_pct'].median():.3f}")


if __name__ == "__main__":
    main()
