"""Day 5: degree-preserving permutation null -> z-scored DWPC -> re-evaluation.

Why: damped path counts still favour compounds that bind many genes (hub audit: top-20
candidates bind ~3-4x more genes than the typical compound). A z-score against a
degree-preserving null asks a stricter question: "more paths to this disease than expected
for a compound/gene/pathway with this degree profile?"

Method (Himmelstein et al. 2017 style, own implementation):
  * XSwap edge swaps preserve every node's degree in each permuted relation
    (CbG, CuG, CdG, DaG, DuG, DdG, GpPW, GiG).
  * For each permutation recompute the four mechanistic metapaths; collect per-pair mean/sd.
  * z = (observed - null mean) / sqrt(null sd^2 + floor^2); floor = 10th percentile of the
    positive null sds (stabilises pairs that are almost always zero in the null).
  * Rank with the coverage-aware percentile aggregation (all metapaths treated as signed).

Note the trade-off it exposes: known indications correlate with compound degree (popular,
well-studied drugs), so removing degree can LOWER benchmark AUROC while making candidate
lists more meaningful. Both are reported.

Outputs: data/processed/null_*.npz (cached), zscores_w<w>.npz;
outputs/tables/eval_z_summary.csv, eval_z_per_disease.csv, ranking_<slug>_z_w<w>.csv,
known_ranks_<slug>_z_w<w>.csv

Usage:
    python3 src/08_permutation_null.py                    # ~20 permutations (a few minutes)
    python3 src/08_permutation_null.py --n-perm 50        # tighter null
    python3 src/08_permutation_null.py --multiplier 5     # faster, less mixing
"""
from __future__ import annotations

import argparse
import time

import numpy as np
import pandas as pd
import scipy.sparse as sp

from dwpc import (MECHANISTIC, SYMMETRIC, aggregate_vector, build_matrices, compute_dwpc,
                  coverage, rank_disease)
from evalutils import metrics, rho, summarize
from utils import DATA_PROC, OUT_TAB, find_nodes, load_processed

PERMUTE = ["CbG", "CuG", "CdG", "DaG", "DuG", "DdG", "GpPW", "GiG"]
BIND = ["CbGaD", "CbGiGaD", "CbGpPWpGaD"]


def permute_relation(M: sp.spmatrix, symmetric: bool, multiplier: float,
                     rng: np.random.Generator) -> sp.csr_matrix:
    """Degree-preserving XSwap on a binary bipartite (or symmetric) adjacency matrix."""
    coo = (sp.triu(M, k=1) if symmetric else M).tocoo()
    edges = list(zip(coo.row.tolist(), coo.col.tolist()))
    m = len(edges)
    S = set(edges)
    n_swaps = int(multiplier * m)
    ii = rng.integers(0, m, n_swaps).tolist()
    jj = rng.integers(0, m, n_swaps).tolist()
    flips = (rng.random(n_swaps) < 0.5).tolist() if symmetric else None
    for t in range(n_swaps):
        i, j = ii[t], jj[t]
        if i == j:
            continue
        a, b = edges[i]
        c, d = edges[j]
        if symmetric and flips[t]:
            c, d = d, c
        if a == c or b == d:
            continue
        if symmetric:
            if a == d or c == b:
                continue
            e1 = (a, d) if a < d else (d, a)
            e2 = (c, b) if c < b else (b, c)
        else:
            e1, e2 = (a, d), (c, b)
        if e1 in S or e2 in S or e1 == e2:
            continue
        S.discard(edges[i])
        S.discard(edges[j])
        edges[i], edges[j] = e1, e2
        S.add(e1)
        S.add(e2)
    rows, cols = zip(*edges)
    P = sp.csr_matrix((np.ones(m), (rows, cols)), shape=M.shape)
    if symmetric:
        P = P + P.T
    P = P.tocsr()
    P.data[:] = 1.0
    return P


def null_stats(A: dict, w: float, n_perm: int, multiplier: float, rng: np.random.Generator):
    s1 = {k: 0.0 for k in MECHANISTIC}
    s2 = {k: 0.0 for k in MECHANISTIC}
    t0 = time.time()
    for p in range(n_perm):
        Ap = dict(A)
        for rel in PERMUTE:
            Ap[rel] = permute_relation(A[rel], rel in SYMMETRIC, multiplier, rng)
        sc = compute_dwpc(Ap, w)
        for k in MECHANISTIC:
            s1[k] = s1[k] + sc[k]
            s2[k] = s2[k] + sc[k] ** 2
        if p == 0:
            print(f"  first permutation took {time.time() - t0:.1f}s "
                  f"-> ~{(time.time() - t0) * n_perm / 60:.1f} min total")
    mu, sd = {}, {}
    for k in MECHANISTIC:
        mu[k] = s1[k] / n_perm
        var = np.maximum(s2[k] / n_perm - mu[k] ** 2, 0.0) * n_perm / max(n_perm - 1, 1)
        sd[k] = np.sqrt(var)
    print(f"  null done in {time.time() - t0:.0f}s")
    return mu, sd


def zscores(obs: dict, mu: dict, sd: dict) -> dict:
    Z = {}
    for k in MECHANISTIC:
        pos = sd[k][sd[k] > 0]
        floor = float(np.percentile(pos, 10)) if len(pos) else 1.0
        Z[k] = (obs[k] - mu[k]) / np.sqrt(sd[k] ** 2 + floor ** 2)
    return Z


def method_vectors(d, obs, Z, cov, deg, min_nonzero):
    out = {"degree_CbG": deg}
    for mp in MECHANISTIC:
        v = obs[mp][:, d].astype(float).copy()
        v[~cov[mp]] = np.nan
        out[f"{mp}@dwpc"] = v
        z = Z[mp][:, d].astype(float).copy()
        z[~cov[mp]] = np.nan
        out[f"{mp}@z"] = z
    for name, mps in (("agg_mech_mean", MECHANISTIC), ("agg_bind_mean", BIND)):
        out[f"{name}@dwpc"] = aggregate_vector(obs, cov, d, mps, "mean", None, min_nonzero)[0]
        out[f"{name}@z"] = aggregate_vector(Z, cov, d, mps, "mean", None, 0, signed_all=True)[0]
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--w", type=float, default=0.6)
    ap.add_argument("--n-perm", type=int, default=20)
    ap.add_argument("--multiplier", type=float, default=10.0)
    ap.add_argument("--min-pos", type=int, default=3)
    ap.add_argument("--min-nonzero", type=int, default=20)
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--recompute", action="store_true")
    ap.add_argument("--diseases", default="alzheimer,parkinson")
    ap.add_argument("--top", type=int, default=25)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    nodes, edges = load_processed()
    mats = build_matrices(nodes, edges)
    A = mats.A
    cov = coverage(A)
    dis = mats.tables["Disease"]
    comp = mats.tables["Compound"]
    deg = np.asarray(A["CbG"].sum(axis=1)).ravel().astype(float)
    CtD = A["CtD"].toarray() > 0
    bench = np.where(CtD.sum(axis=0) >= args.min_pos)[0]
    obs = compute_dwpc(A, args.w)

    cache = DATA_PROC / f"null_w{args.w:g}_p{args.n_perm}_m{args.multiplier:g}.npz"
    if cache.exists() and not args.recompute:
        f = np.load(cache)
        mu = {k: f[f"mu_{k}"] for k in MECHANISTIC}
        sd = {k: f[f"sd_{k}"] for k in MECHANISTIC}
        print(f"Loaded cached null: {cache.name}")
    else:
        print(f"Building null: {args.n_perm} permutations, multiplier {args.multiplier:g}, w={args.w:g}")
        mu, sd = null_stats(A, args.w, args.n_perm, args.multiplier, rng)
        np.savez_compressed(cache, **{f"mu_{k}": mu[k] for k in MECHANISTIC},
                            **{f"sd_{k}": sd[k] for k in MECHANISTIC})
    Z = zscores(obs, mu, sd)
    np.savez_compressed(DATA_PROC / f"zscores_w{args.w:g}.npz", **Z,
                        compound_ids=comp["id"].to_numpy(), disease_ids=dis["id"].to_numpy())

    # ---------------- pan-disease benchmark ----------------
    rows = []
    for d in bench:
        pos = CtD[:, d]
        for name, v in method_vectors(d, obs, Z, cov, deg, args.min_nonzero).items():
            auc, med, h5, h10 = metrics(v, pos)
            rows.append(dict(disease_id=dis.loc[d, "id"], disease=dis.loc[d, "name"], n_pos=int(pos.sum()),
                             method=name, auroc=auc, median_rank_pct=med, hit5=h5, hit10=h10,
                             rho_degree=rho(v, deg)))
    per = pd.DataFrame(rows)
    per.to_csv(OUT_TAB / "eval_z_per_disease.csv", index=False)
    summ = summarize(per, args.boot, rng)
    summ.to_csv(OUT_TAB / "eval_z_summary.csv", index=False)
    print(f"\n=== Pan-disease benchmark, w={args.w:g}: {len(bench)} diseases, {int(CtD[:, bench].sum())} positives ===")
    print(summ[["method", "mean_auroc", "ci_lo", "ci_hi", "wins_vs_degree", "p_vs_degree",
                "mean_hit5", "rho_degree"]].round(3).to_string(index=False))

    # ---------------- AD / PD rankings ----------------
    for query in [q.strip() for q in args.diseases.split(",") if q.strip()]:
        hits = find_nodes(nodes, query, kind="Disease")
        if hits.empty:
            continue
        row = hits.iloc[0]
        d = mats.local_index(int(row["idx"]))
        full, used = rank_disease(mats, Z, cov, d, MECHANISTIC, mode="mean",
                                  min_nonzero=0, signed_all=True)
        full = full.rename(columns={f"dwpc_{mp}": f"z_{mp}" for mp in MECHANISTIC})
        ci = full["cidx"].to_numpy()
        for mp in used:
            full[f"obs_{mp}"] = obs[mp][ci, d]
            full[f"null_mu_{mp}"] = mu[mp][ci, d]

        print(f"\n{'=' * 78}\n{row['name']} - z-scored ranking (w={args.w:g}, {args.n_perm} permutations)\n{'=' * 78}")
        cols = ["rank", "compound", "score", "n_covered"] + [f"pct_{m}" for m in used] + [f"z_{m}" for m in used]
        cand = full[~full["known"] & full["score"].notna()]
        show = cand.head(args.top)[cols].copy()
        show["rank"] = show["rank"].astype(int)
        print(show.round(2).to_string(index=False))
        top = cand.head(20)
        print(f"\nHub audit: median # bound genes - top-20 {np.median(deg[top['cidx']]):.0f} "
              f"vs all covered compounds {np.median(deg[cov['CbGaD']]):.0f}")

        ranked = full[full["score"].notna()]
        known = full[full["known"] & full["score"].notna()].copy()
        known["n_ranked"] = len(ranked)
        known["rank_pct"] = known["rank"] / len(ranked)
        print(f"Median rank percentile of known drugs: {known['rank_pct'].median():.3f} "
              f"(n={len(known)}, random ~0.5)")

        slug = query.lower().replace(" ", "_")
        cand.head(200).to_csv(OUT_TAB / f"ranking_{slug}_z_w{args.w:g}.csv", index=False)
        known[["compound", "rank", "n_ranked", "rank_pct", "score"] + [f"z_{m}" for m in used]].to_csv(
            OUT_TAB / f"known_ranks_{slug}_z_w{args.w:g}.csv", index=False)


if __name__ == "__main__":
    main()
