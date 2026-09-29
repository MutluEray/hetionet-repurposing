"""Day 5b: consensus candidates - DWPC (signal) AND z-score (specificity).

Findings that motivate this script (pan-disease benchmark, w=0.6):
  * DWPC on the binding metapaths is the best scorer (AUROC ~0.79) but favours compounds
    that bind many genes.
  * The permutation z-score removes that degree effect (Spearman with degree < 0) at the
    cost of AUROC (~0.70, about the degree baseline) - so a degree-independent signal
    exists, but it is not free.
  * SigNet (expression-signature reversal) does not improve the benchmark; it is kept as
    an ANNOTATION (it is the only evidence with a direction of effect), not a scored term.

Consensus = mean of the two rank percentiles (DWPC-binding aggregate, z-binding aggregate);
"consensus_min" = the lower of the two (must be high on both). Both are benchmarked against
the degree baseline. Candidate tables include a blank `triage` column for the manual
literature step (already trialled / risk factor / endogenous / diagnostic / cytotoxic ...).

Run 08_permutation_null.py first (it writes data/processed/zscores_w<w>.npz).

Outputs: outputs/tables/eval_consensus_summary.csv, candidates_<slug>.csv

Usage:
    python3 src/09_consensus_candidates.py
    python3 src/09_consensus_candidates.py --top 30
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from dwpc import MECHANISTIC, aggregate_vector, build_matrices, compute_dwpc, coverage, known_mask
from evalutils import metrics, rho, summarize
from utils import DATA_PROC, OUT_TAB, find_nodes, load_processed

BIND = ["CbGaD", "CbGiGaD", "CbGpPWpGaD"]


def pct(v: np.ndarray) -> np.ndarray:
    out = np.full(len(v), np.nan)
    ok = np.isfinite(v)
    out[ok] = rankdata(v[ok], method="average") / ok.sum()
    return out


def vectors(obs, Z, cov, d, min_nonzero):
    s_d = aggregate_vector(obs, cov, d, BIND, "mean", None, min_nonzero)[0]
    s_z = aggregate_vector(Z, cov, d, BIND, "mean", None, 0, signed_all=True)[0]
    p_d, p_z = pct(s_d), pct(s_z)
    return {"dwpc_bind": s_d, "z_bind": s_z, "consensus_mean": (p_d + p_z) / 2,
            "consensus_min": np.minimum(p_d, p_z)}, p_d, p_z


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--w", type=float, default=0.6)
    ap.add_argument("--min-pos", type=int, default=3)
    ap.add_argument("--min-nonzero", type=int, default=20)
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--top", type=int, default=25)
    ap.add_argument("--keep", type=int, default=100, help="rows saved per disease")
    ap.add_argument("--diseases", default="alzheimer,parkinson")
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    zpath = DATA_PROC / f"zscores_w{args.w:g}.npz"
    if not zpath.exists():
        raise SystemExit(f"{zpath} not found - run 08_permutation_null.py first.")
    f = np.load(zpath)
    Z = {k: f[k] for k in MECHANISTIC}

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

    # ---------------- benchmark ----------------
    rows = []
    for d in bench:
        pos = CtD[:, d]
        vecs, _, _ = vectors(obs, Z, cov, d, args.min_nonzero)
        vecs["degree_CbG"] = deg
        for name, v in vecs.items():
            auc, med, h5, h10 = metrics(v, pos)
            rows.append(dict(disease_id=dis.loc[d, "id"], disease=dis.loc[d, "name"], n_pos=int(pos.sum()),
                             method=name, auroc=auc, median_rank_pct=med, hit5=h5, hit10=h10,
                             rho_degree=rho(v, deg)))
    per = pd.DataFrame(rows)
    summ = summarize(per, args.boot, rng)
    summ.to_csv(OUT_TAB / "eval_consensus_summary.csv", index=False)
    print(f"=== Pan-disease benchmark, w={args.w:g}: {len(bench)} diseases ===")
    print(summ[["method", "mean_auroc", "ci_lo", "ci_hi", "wins_vs_degree", "p_vs_degree",
                "mean_hit5", "rho_degree"]].round(3).to_string(index=False))

    # ---------------- AD / PD consensus candidates ----------------
    for query in [q.strip() for q in args.diseases.split(",") if q.strip()]:
        hits = find_nodes(nodes, query, kind="Disease")
        if hits.empty:
            continue
        row = hits.iloc[0]
        d = mats.local_index(int(row["idx"]))
        vecs, p_d, p_z = vectors(obs, Z, cov, d, args.min_nonzero)
        cons = vecs["consensus_mean"]
        sig_cov = cov["SigNet"]
        df = pd.DataFrame({
            "compound_id": comp["id"], "compound": comp["name"], "cidx": np.arange(len(comp)),
            "consensus": cons, "pct_dwpc": p_d, "pct_z": p_z,
            "n_bound_genes": deg,
            "z_CbGaD": Z["CbGaD"][:, d], "z_CbGiGaD": Z["CbGiGaD"][:, d], "z_CbGpPWpGaD": Z["CbGpPWpGaD"][:, d],
            "obs_CbGaD": obs["CbGaD"][:, d],
            "sig_reversal": np.where(sig_cov, obs["SigRev"][:, d], np.nan),
            "sig_concordance": np.where(sig_cov, obs["SigConc"][:, d], np.nan),
            "known": known_mask(A, d), "triage": "",
        })
        df["rank"] = df["consensus"].rank(ascending=False, method="min")
        df = df.sort_values("consensus", ascending=False, na_position="last").reset_index(drop=True)

        cand = df[~df["known"] & df["consensus"].notna()]
        print(f"\n{'=' * 78}\n{row['name']} - consensus candidates (known drugs excluded)\n{'=' * 78}")
        cols = ["rank", "compound", "consensus", "pct_dwpc", "pct_z", "n_bound_genes",
                "z_CbGaD", "z_CbGpPWpGaD", "sig_reversal", "sig_concordance"]
        show = cand.head(args.top)[cols].copy()
        show["rank"] = show["rank"].astype(int)
        print(show.round(2).to_string(index=False))
        top = cand.head(20)
        print(f"\nHub audit: median # bound genes - top-20 {top['n_bound_genes'].median():.0f} "
              f"vs all covered compounds {np.median(deg[cov['CbGaD']]):.0f}")
        ranked = df[df["consensus"].notna()]
        known = df[df["known"] & df["consensus"].notna()]
        print(f"Known drugs: median rank percentile {(known['rank'] / len(ranked)).median():.3f} "
              f"(n={len(known)}, random ~0.5)")
        slug = query.lower().replace(" ", "_")
        cand.head(args.keep).to_csv(OUT_TAB / f"candidates_{slug}.csv", index=False)


if __name__ == "__main__":
    main()
