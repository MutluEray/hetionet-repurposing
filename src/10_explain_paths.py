"""Day 5c: explanation layer - explicit paths, contribution shares and counterfactual ablation.

For the shortlisted candidates of each disease (top consensus candidates after removing
triage categories such as 'endogenous', 'cytotoxic_or_metal', 'risk_factor', ...), this script:

  1. retrieves every supporting path with SPARQL (rdflib, core RDF graph from 04_build_rdf.py):
        CbGaD       compound -binds-> gene <-associates- disease
        CbGiGaD     compound -binds-> gene -interacts- gene <-associates- disease
        CbGpPWpGaD  compound -binds-> gene -participates-> pathway <-participates- gene <-associates- disease
  2. weights each path exactly as the matrix engine does (product of (deg_src*deg_dst)^-w per edge)
     and checks that the path weights sum to the matrix DWPC (built-in consistency check),
  3. reports which disease genes / pathways drive each candidate (contribution shares),
  4. runs a counterfactual: remove the candidate's top-1 / top-3 driving disease genes from the
     disease's gene set and recompute the DWPC-binding rank (z-scores are not recomputed).

Inputs : outputs/tables/candidates_<slug>.csv (09), outputs/rdf/hetionet_core.nt (04),
         data/triage.csv (optional: compound, category, note)
Outputs: outputs/tables/paths_<slug>.csv, explain_summary_<slug>.csv, candidates_<slug>_triaged.csv

Usage:
    python3 src/10_explain_paths.py
    python3 src/10_explain_paths.py --k 8 --exclude-categories endogenous,diagnostic
"""
from __future__ import annotations

import argparse
import time

import numpy as np
import pandas as pd
from rdflib import Graph

from dwpc import aggregate_vector, build_matrices, compute_dwpc, coverage
from utils import OUT_RDF, OUT_TAB, ROOT, find_nodes, load_processed, node_uri

BIND = ["CbGaD", "CbGiGaD", "CbGpPWpGaD"]
DEFAULT_EXCLUDE = "endogenous,diagnostic,cytotoxic_or_metal,risk_factor,wrong_direction_likely"

PFX = "PREFIX rel: <http://het.io/rel/>\n"
Q_A = "SELECT DISTINCT ?g WHERE {{ {C} rel:CbG ?g . {D} rel:DaG ?g }}"
Q_B = ("SELECT DISTINCT ?g1 ?g2 WHERE {{ {C} rel:CbG ?g1 . ?g1 (rel:GiG|^rel:GiG) ?g2 . "
       "{D} rel:DaG ?g2 . FILTER(?g1 != ?g2) }}")
Q_C = ("SELECT DISTINCT ?g1 ?p ?g2 WHERE {{ {C} rel:CbG ?g1 . ?g1 rel:GpPW ?p . ?g2 rel:GpPW ?p . "
       "{D} rel:DaG ?g2 . FILTER(?g1 != ?g2) }}")


def relation_degrees(A: dict) -> dict:
    return {r: (np.asarray(A[r].sum(axis=1)).ravel(), np.asarray(A[r].sum(axis=0)).ravel())
            for r in ("CbG", "DaG", "GpPW", "GiG")}


def make_ew(deg: dict, w: float):
    def ew(rel: str, u: int, v: int) -> float:
        dr, dc = deg[rel]
        return float((dr[u] * dc[v]) ** (-w))
    return ew


def bind_rank(scores: dict, cov: dict, d: int, c: int, min_nonzero: int) -> tuple[int, int]:
    s, _, _ = aggregate_vector(scores, cov, d, BIND, "mean", None, min_nonzero)
    ok = np.isfinite(s)
    if not ok[c]:
        return -1, int(ok.sum())
    return int(1 + (s[ok] > s[c]).sum()), int(ok.sum())


def ablate(A: dict, d: int, genes: list[int]) -> dict:
    M = A["DaG"].tolil()
    for g in genes:
        M[d, g] = 0
    M = M.tocsr()
    M.eliminate_zeros()
    A2 = dict(A)
    A2["DaG"] = M
    return A2


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--w", type=float, default=0.6)
    ap.add_argument("--k", type=int, default=8, help="shortlisted candidates per disease")
    ap.add_argument("--min-nonzero", type=int, default=20)
    ap.add_argument("--diseases", default="alzheimer,parkinson")
    ap.add_argument("--exclude-categories", default=DEFAULT_EXCLUDE)
    ap.add_argument("--graph", default=str(OUT_RDF / "hetionet_core.nt"))
    ap.add_argument("--max-paths", type=int, default=50, help="paths saved per candidate and metapath")
    args = ap.parse_args()
    excl = {x.strip() for x in args.exclude_categories.split(",") if x.strip()}

    nodes, edges = load_processed()
    mats = build_matrices(nodes, edges)
    A = mats.A
    cov = coverage(A)
    comp, gene, path, dis = (mats.tables[k] for k in ("Compound", "Gene", "Pathway", "Disease"))
    obs = compute_dwpc(A, args.w)
    deg = relation_degrees(A)
    ew = make_ew(deg, args.w)
    n_bound = np.asarray(A["CbG"].sum(axis=1)).ravel()

    t = time.time()
    g = Graph()
    g.parse(args.graph, format="nt")
    print(f"Loaded {len(g):,} triples in {time.time() - t:.1f}s")
    gene_loc = {node_uri(i): k for k, i in enumerate(gene["id"])}
    path_loc = {node_uri(i): k for k, i in enumerate(path["id"])}

    triage_path = ROOT / "data" / "triage.csv"
    triage = (pd.read_csv(triage_path).assign(key=lambda x: x["compound"].str.lower())
              .set_index("key")[["category", "note"]] if triage_path.exists() else None)

    for query in [q.strip() for q in args.diseases.split(",") if q.strip()]:
        hits = find_nodes(nodes, query, kind="Disease")
        cpath = OUT_TAB / f"candidates_{query.lower()}.csv"
        if hits.empty or not cpath.exists():
            print(f"\nSkipping '{query}' (disease or {cpath.name} missing - run 09 first)")
            continue
        row = hits.iloc[0]
        d = mats.local_index(int(row["idx"]))
        D = f"<{node_uri(row['id'])}>"
        cand = pd.read_csv(cpath)
        cand["category"], cand["note"] = "", ""
        if triage is not None:
            key = cand["compound"].str.lower()
            cand["category"] = key.map(triage["category"]).fillna("")
            cand["note"] = key.map(triage["note"]).fillna("")
        cand["shortlisted"] = False
        short_idx = cand.index[~cand["category"].isin(excl)][: args.k]
        cand.loc[short_idx, "shortlisted"] = True
        slug = query.lower()
        cand.head(50).to_csv(OUT_TAB / f"candidates_{slug}_triaged.csv", index=False)
        print(f"\n{'=' * 78}\n{row['name']}: explaining {len(short_idx)} shortlisted candidates "
              f"(excluded categories: {sorted(excl)})\n{'=' * 78}")
        base_scores = obs
        path_rows, sum_rows = [], []

        for i in short_idx:
            r = cand.loc[i]
            c = int(r["cidx"])
            C = f"<{node_uri(r['compound_id'])}>"
            cname = r["compound"]
            paths = {"CbGaD": [], "CbGiGaD": [], "CbGpPWpGaD": []}
            for (gg,) in g.query(PFX + Q_A.format(C=C, D=D)):
                gi = gene_loc[str(gg)]
                paths["CbGaD"].append(((gi,), ew("CbG", c, gi) * ew("DaG", d, gi), gi))
            for g1, g2 in g.query(PFX + Q_B.format(C=C, D=D)):
                a, b = gene_loc[str(g1)], gene_loc[str(g2)]
                paths["CbGiGaD"].append(((a, b), ew("CbG", c, a) * ew("GiG", a, b) * ew("DaG", d, b), b))
            for g1, p, g2 in g.query(PFX + Q_C.format(C=C, D=D)):
                a, pw, b = gene_loc[str(g1)], path_loc[str(p)], gene_loc[str(g2)]
                paths["CbGpPWpGaD"].append(
                    ((a, pw, b), ew("CbG", c, a) * ew("GpPW", a, pw) * ew("GpPW", b, pw) * ew("DaG", d, b), b))

            totals, check_ok = {}, True
            share_gene: dict[int, float] = {}
            share_path: dict[int, float] = {}
            for mp, plist in paths.items():
                tot = sum(w_ for _, w_, _ in plist)
                totals[mp] = tot
                ref = float(obs[mp][c, d])
                if not np.isclose(tot, ref, rtol=1e-6, atol=1e-9):
                    check_ok = False
                if tot > 0:
                    per_gene: dict[int, float] = {}
                    for nodes_t, w_, end in plist:
                        per_gene[end] = per_gene.get(end, 0.0) + w_
                        if mp == "CbGpPWpGaD":
                            share_path[nodes_t[1]] = share_path.get(nodes_t[1], 0.0) + w_ / tot
                    for k_, v_ in per_gene.items():
                        share_gene[k_] = share_gene.get(k_, 0.0) + v_ / tot / 3.0
                for nodes_t, w_, end in sorted(plist, key=lambda x: -x[1])[: args.max_paths]:
                    if mp == "CbGaD":
                        desc = f"{cname} -binds-> {gene.loc[nodes_t[0], 'name']} <-assoc- {row['name']}"
                    elif mp == "CbGiGaD":
                        desc = (f"{cname} -binds-> {gene.loc[nodes_t[0], 'name']} -interacts- "
                                f"{gene.loc[nodes_t[1], 'name']} <-assoc- {row['name']}")
                    else:
                        desc = (f"{cname} -binds-> {gene.loc[nodes_t[0], 'name']} -in-> "
                                f"{path.loc[nodes_t[1], 'name']} <-in- {gene.loc[nodes_t[2], 'name']} "
                                f"<-assoc- {row['name']}")
                    path_rows.append(dict(disease=row["name"], compound=cname, metapath=mp,
                                          weight=w_, share_in_metapath=w_ / tot if tot else np.nan, path=desc))

            top_genes = sorted(share_gene, key=share_gene.get, reverse=True)[:3]
            top_pw = max(share_path, key=share_path.get) if share_path else None
            rank0, n_ranked = bind_rank(obs, cov, d, c, args.min_nonzero)
            r1 = r3 = np.nan
            if top_genes:
                r1 = bind_rank(compute_dwpc(ablate(A, d, top_genes[:1]), args.w), cov, d, c, args.min_nonzero)[0]
                r3 = bind_rank(compute_dwpc(ablate(A, d, top_genes[:3]), args.w), cov, d, c, args.min_nonzero)[0]
            sum_rows.append(dict(
                compound=cname, category=r["category"], consensus=r["consensus"],
                n_bound_genes=int(n_bound[c]),
                **{f"n_paths_{mp}": len(paths[mp]) for mp in BIND},
                **{f"dwpc_{mp}": totals[mp] for mp in BIND}, sum_check_ok=check_ok,
                top_gene=gene.loc[top_genes[0], "name"] if top_genes else "",
                top_gene_share=share_gene[top_genes[0]] if top_genes else np.nan,
                top3_genes=", ".join(gene.loc[x, "name"] for x in top_genes),
                top_pathway=path.loc[top_pw, "name"] if top_pw is not None else "",
                top_pathway_share=share_path[top_pw] if top_pw is not None else np.nan,
                dwpc_rank=rank0, rank_no_top1=r1, rank_no_top3=r3, n_ranked=n_ranked))

        summ = pd.DataFrame(sum_rows)
        pd.DataFrame(path_rows).to_csv(OUT_TAB / f"paths_{slug}.csv", index=False)
        summ.to_csv(OUT_TAB / f"explain_summary_{slug}.csv", index=False)
        show = summ[["compound", "category", "n_bound_genes", "top_gene", "top_gene_share", "top_pathway",
                     "top_pathway_share", "dwpc_rank", "rank_no_top1", "rank_no_top3", "sum_check_ok"]]
        print(show.round(2).to_string(index=False))
        top_all = pd.Series([x for s in summ["top3_genes"] for x in s.split(", ") if x]).value_counts()
        print(f"\nMost frequent driving disease genes across the shortlist:\n{top_all.head(8).to_string()}")
        pr = pd.DataFrame(path_rows)
        for cname in summ["compound"].head(3):
            sub = pr[(pr["compound"] == cname)].sort_values("weight", ascending=False).head(4)
            print(f"\nTop paths for {cname}:")
            for _, x in sub.iterrows():
                print(f"  [{x['metapath']}, share {x['share_in_metapath']:.2f}] {x['path']}")


if __name__ == "__main__":
    main()
