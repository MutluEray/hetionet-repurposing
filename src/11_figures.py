"""Day 6: portfolio figures from the pipeline outputs (v2).

Reads the tables written by 07-10 and draws:
    fig_benchmark.png         AUROC (dot + 95% CI) and the AUROC-vs-degree-correlation trade-off
    fig_damping.png           damping sweep: AUROC and degree correlation vs w
    fig_hub_audit.png         # bound genes: all compounds vs top-20 consensus candidates
    fig_ablation.png          rank before/after removing the top driving disease gene
    evidence_<slug>.png       evidence card table for the shortlisted candidates
    paths_<slug>_panel.png    path graphs of the top-4 shortlisted candidates
    path_<slug>_<compound>.png  one path graph per shortlisted candidate

Support tier (from the gene ablation): 'multi-gene' if the candidate stays within the top 5% of
DWPC-binding ranks after its top driving disease gene is removed, otherwise 'single-gene-dependent'.
This measures dependence on one gene, not overall evidence breadth (see the 'w/o top 3' column).

Pathway paths that run through pathways named after a disease (e.g. 'Alzheimers Disease') are hidden
from the path graphs: they are near-tautological as explanations (07b shows they do not change the ranking).

Path graphs: genes that are both bound and disease-associated appear in both roles; nodes are
ordered by a barycenter heuristic to reduce edge crossings.

Usage:
    python3 src/11_figures.py
"""
from __future__ import annotations

import re
import textwrap

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from dwpc import build_matrices
from utils import OUT_FIG, OUT_TAB, load_processed

PARSE = re.compile(r" -binds-> | <-assoc- | -interacts- | -in-> | <-in- ")
DISEASES = ["alzheimer", "parkinson"]
FAMILY_COLOR = {"baseline": "#7f7f7f", "dwpc": "#1f77b4", "z": "#ff7f0e", "consensus": "#2ca02c"}
NODE_COLOR = {"compound": "#4c78a8", "bgene": "#9ecae9", "pathway": "#a1d99b",
              "dgene": "#fdae6b", "disease": "#e6550d"}
EDGE_COLOR = {"binds": "#1f77b4", "interacts": "#7f7f7f", "in": "#2ca02c", "assoc": "#d62728"}
EDGE_LABEL = {"binds": "binds", "interacts": "interacts with", "in": "participates in", "assoc": "associated with"}
LAYER = {"compound": 0, "bgene": 1, "pathway": 2, "dgene": 3, "disease": 4}
TIER_OK, TIER_BAD = "multi-gene", "single-gene-dependent"


def read(name: str) -> pd.DataFrame | None:
    p = OUT_TAB / name
    if not p.exists():
        print(f"  (missing {name} - skipping dependent figure)")
        return None
    return pd.read_csv(p)


# --------------------------------------------------------------------------- benchmark
def fig_benchmark() -> None:
    z = read("eval_z_summary.csv")
    cons = read("eval_consensus_summary.csv")
    if z is None:
        return
    df = z.copy()
    if cons is not None:
        df = pd.concat([df, cons[cons["method"].isin(["consensus_mean", "consensus_min"])]], ignore_index=True)
    names = {
        "degree_CbG": ("Degree baseline (# bound genes)", "baseline"),
        "agg_bind_mean@dwpc": ("DWPC, binding metapaths", "dwpc"),
        "agg_mech_mean@dwpc": ("DWPC, + signature", "dwpc"),
        "agg_bind_mean@z": ("z-score, binding metapaths", "z"),
        "agg_mech_mean@z": ("z-score, + signature", "z"),
        "consensus_mean": ("Consensus (mean of DWPC & z pct)", "consensus"),
        "consensus_min": ("Consensus (min of DWPC & z pct)", "consensus"),
        "CbGaD@dwpc": ("CbGaD (DWPC)", "dwpc"), "CbGiGaD@dwpc": ("CbGiGaD (DWPC)", "dwpc"),
        "CbGpPWpGaD@dwpc": ("CbGpPWpGaD (DWPC)", "dwpc"), "SigNet@dwpc": ("SigNet (DWPC)", "dwpc"),
        "CbGaD@z": ("CbGaD (z)", "z"), "CbGiGaD@z": ("CbGiGaD (z)", "z"),
        "CbGpPWpGaD@z": ("CbGpPWpGaD (z)", "z"), "SigNet@z": ("SigNet (z)", "z"),
    }
    df = df[df["method"].isin(names)].sort_values("mean_auroc")
    fig, (a, b) = plt.subplots(1, 2, figsize=(15.5, 6.5), gridspec_kw={"width_ratios": [1.1, 1]})
    y = np.arange(len(df))
    for yi, (_, r) in zip(y, df.iterrows()):
        c = FAMILY_COLOR[names[r["method"]][1]]
        a.plot([r["ci_lo"], r["ci_hi"]], [yi, yi], color=c, lw=2.2)
        a.plot(r["mean_auroc"], yi, "o", color=c, ms=8, mec="black", mew=0.6)
    a.set_yticks(y)
    a.set_yticklabels([names[m][0] for m in df["method"]], fontsize=8)
    if (df["method"] == "degree_CbG").any():
        a.axvline(float(df.loc[df["method"] == "degree_CbG", "mean_auroc"].iloc[0]), color="gray", ls="--", lw=1)
    a.axvline(0.5, color="black", lw=0.6)
    a.set_xlim(0.5, 0.86)
    a.grid(axis="x", alpha=0.25)
    a.set_xlabel("Macro AUROC across diseases (known 'treats' edges; 0.5 = chance)")
    a.set_title("Pan-disease benchmark (dot = mean, line = 95% bootstrap CI over diseases)")
    a.legend(handles=[Line2D([0], [0], marker="o", color=c, lw=2, label=l, mec="black")
                      for l, c in (("degree baseline", FAMILY_COLOR["baseline"]), ("DWPC", FAMILY_COLOR["dwpc"]),
                                   ("z-score vs permutation null", FAMILY_COLOR["z"]),
                                   ("consensus", FAMILY_COLOR["consensus"]))],
             loc="lower right", fontsize=8)

    agg = df[df["method"].isin(["degree_CbG", "agg_bind_mean@dwpc", "agg_mech_mean@dwpc", "agg_bind_mean@z",
                                "agg_mech_mean@z", "consensus_mean", "consensus_min"])]
    offsets = {"degree_CbG": (-150, 12), "agg_bind_mean@dwpc": (10, 10), "agg_mech_mean@dwpc": (10, -16),
               "agg_bind_mean@z": (10, 8), "agg_mech_mean@z": (10, -14),
               "consensus_mean": (-10, 8), "consensus_min": (-10, -14)}
    left = {"consensus_mean", "consensus_min"}
    for _, r in agg.iterrows():
        fam = names[r["method"]][1]
        b.errorbar(r["rho_degree"], r["mean_auroc"], yerr=[[r["mean_auroc"] - r["ci_lo"]], [r["ci_hi"] - r["mean_auroc"]]],
                   fmt="o", color=FAMILY_COLOR[fam], capsize=3, ms=8, mec="black", mew=0.6)
        b.annotate(names[r["method"]][0], (r["rho_degree"], r["mean_auroc"]), textcoords="offset points",
                   xytext=offsets.get(r["method"], (8, 6)), fontsize=8,
                   ha="right" if r["method"] in left else "left")
    b.axhline(0.5, color="black", lw=0.6)
    b.axvline(0, color="gray", ls=":", lw=1)
    b.set_xlim(-0.75, 1.25)
    b.grid(alpha=0.25)
    b.set_xlabel("Mean Spearman correlation of score with compound degree (# bound genes)")
    b.set_ylabel("Macro AUROC")
    b.set_title("Performance vs dependence on # bound genes")
    fig.tight_layout()
    fig.savefig(OUT_FIG / "fig_benchmark.png", dpi=200)
    plt.close(fig)
    print("  fig_benchmark.png")


def fig_damping() -> None:
    d = read("eval_damping.csv")
    if d is None:
        return
    fig, (a, b) = plt.subplots(1, 2, figsize=(13, 4.8))
    for m, g in d.groupby("method"):
        a.plot(g["w"], g["mean_auroc"], marker="o", label=m)
        b.plot(g["w"], g["mean_rho_degree"], marker="o", label=m)
    a.set_xlabel("damping exponent w")
    a.set_ylabel("Macro AUROC")
    a.set_title("Accuracy is flat in w ...")
    b.set_xlabel("damping exponent w")
    b.set_ylabel("Spearman(score, # bound genes)")
    b.set_title("... but dependence on degree falls")
    a.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_FIG / "fig_damping.png", dpi=200)
    plt.close(fig)
    print("  fig_damping.png")


def fig_hub_audit() -> None:
    nodes, edges = load_processed()
    A = build_matrices(nodes, edges).A
    deg = np.asarray(A["CbG"].sum(axis=1)).ravel()
    data, labels = [deg[deg > 0]], ["All compounds\nwith binding data"]
    for q in DISEASES:
        c = read(f"candidates_{q}.csv")
        if c is not None:
            data.append(c.head(20)["n_bound_genes"].to_numpy())
            labels.append(f"Top-20 consensus\n{q.capitalize()}'s")
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.boxplot(data, showfliers=True)
    ax.set_xticklabels(labels)
    ax.set_yscale("log")
    ax.set_ylabel("# genes bound by the compound")
    ax.set_title("Hub audit")
    fig.tight_layout()
    fig.savefig(OUT_FIG / "fig_hub_audit.png", dpi=200)
    plt.close(fig)
    print("  fig_hub_audit.png")


# --------------------------------------------------------------------------- candidates
def tier(summ: pd.DataFrame) -> pd.DataFrame:
    s = summ.copy()
    s["tier"] = np.where(s["rank_no_top1"] <= 0.05 * s["n_ranked"], TIER_OK, TIER_BAD)
    return s


def fig_ablation(all_summ: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(7.8, 6.2))
    for disease, g in all_summ.groupby("disease"):
        for _, r in g.iterrows():
            col = "#2ca02c" if r["tier"] == TIER_OK else "#d62728"
            ax.scatter(r["dwpc_rank"], r["rank_no_top1"], color=col, marker="o" if disease == "alzheimer" else "s",
                       s=50, edgecolors="black", linewidths=0.5)
            ax.annotate(r["compound"], (r["dwpc_rank"], r["rank_no_top1"]), textcoords="offset points",
                        xytext=(4, 3), fontsize=7)
    m = all_summ["n_ranked"].max()
    ax.plot([1, m], [1, m], color="gray", lw=0.8)
    ax.axhline(0.05 * m, color="gray", ls="--", lw=0.8)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("DWPC-binding rank (observed)")
    ax.set_ylabel("rank after removing the top driving disease gene")
    ax.set_title("Counterfactual gene ablation (circle = AD, square = PD)\n"
                 f"green = {TIER_OK} (stays in top 5%), red = {TIER_BAD}")
    fig.tight_layout()
    fig.savefig(OUT_FIG / "fig_ablation.png", dpi=200)
    plt.close(fig)
    print("  fig_ablation.png")


def fig_evidence(slug: str, summ: pd.DataFrame, cand: pd.DataFrame) -> None:
    m = summ.merge(cand[["compound", "pct_dwpc", "pct_z"]], on="compound", how="left")
    cols = ["compound", "consensus", "pct_dwpc", "pct_z", "n_bound_genes", "top_gene", "top_gene_share",
            "top_pathway", "dwpc_rank", "rank_no_top1", "rank_no_top3", "tier"]
    t = m[cols].copy()
    for c in ("consensus", "pct_dwpc", "pct_z", "top_gene_share"):
        t[c] = t[c].map(lambda x: f"{x:.2f}")
    t["top_pathway"] = t["top_pathway"].map(lambda s: textwrap.shorten(str(s), 34, placeholder="..."))
    for c in ("rank_no_top1", "rank_no_top3"):
        t[c] = t[c].map(lambda x: f"{int(x)}")
    fig, ax = plt.subplots(figsize=(17, 0.55 * len(t) + 1.6))
    ax.axis("off")
    tab = ax.table(cellText=t.values, colLabels=["compound", "consensus", "pct DWPC", "pct z", "# bound genes",
                                                 "top gene", "gene share", "top pathway", "DWPC rank",
                                                 "rank w/o top gene", "rank w/o top 3", "support"],
                   loc="center", cellLoc="left")
    tab.auto_set_font_size(False)
    tab.set_fontsize(8)
    tab.scale(1, 1.5)
    for i, tr in enumerate(t["tier"], start=1):
        tab[(i, len(cols) - 1)].set_facecolor("#c7e9c0" if tr == TIER_OK else "#fcbba1")
    ax.set_title(f"Evidence card - {slug.capitalize()}'s disease shortlist", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT_FIG / f"evidence_{slug}.png", dpi=200)
    plt.close(fig)
    print(f"  evidence_{slug}.png")


def build_graph(sub: pd.DataFrame):
    G = nx.DiGraph()
    for _, r in sub.iterrows():
        parts = PARSE.split(r["path"])
        w = float(r["share_in_metapath"])
        if r["metapath"] == "CbGaD":
            c, g, d = parts
            C, Gn, D = (c, "compound"), (g, "dgene"), (d, "disease")
            chain = [(C, Gn, "binds"), (Gn, D, "assoc")]
        elif r["metapath"] == "CbGiGaD":
            c, g1, g2, d = parts
            C, G1, G2, D = (c, "compound"), (g1, "bgene"), (g2, "dgene"), (d, "disease")
            chain = [(C, G1, "binds"), (G1, G2, "interacts"), (G2, D, "assoc")]
        else:
            c, g1, p, g2, d = parts
            C, G1, P, G2, D = (c, "compound"), (g1, "bgene"), (p, "pathway"), (g2, "dgene"), (d, "disease")
            chain = [(C, G1, "binds"), (G1, P, "in"), (P, G2, "in"), (G2, D, "assoc")]
        for u, v, rel in chain:
            if G.has_edge(u, v):
                G.edges[u, v]["w"] += w
            else:
                G.add_edge(u, v, rel=rel, w=w)
    return G


def layout(G):
    layers = {l: [] for l in range(5)}
    for n in G.nodes:
        layers[LAYER[n[1]]].append(n)
    y = {}

    def assign():
        for ns in layers.values():
            k = len(ns)
            for i, n in enumerate(ns):
                y[n] = (k - 1) / 2 - i

    assign()
    for _ in range(4):
        for l in range(1, 5):
            bary = {n: (np.mean([y[p] for p in G.predecessors(n)]) if G.in_degree(n) else y[n]) for n in layers[l]}
            layers[l].sort(key=lambda n: -bary[n])
            assign()
        for l in range(3, -1, -1):
            bary = {n: (np.mean([y[s] for s in G.successors(n)]) if G.out_degree(n) else y[n]) for n in layers[l]}
            layers[l].sort(key=lambda n: -bary[n])
            assign()
    return {n: (LAYER[n[1]], y[n]) for n in G.nodes}


def rad_for(span: int, y_src: float, y_dst: float) -> float:
    """Curvature of an edge: long (skip-layer) edges bow away from the middle of the drawing."""
    if span >= 3:
        return 0.3 if y_dst < y_src else -0.3
    return 0.05 + 0.07 * span


def draw_candidate(ax, sub: pd.DataFrame, title: str) -> None:
    G = build_graph(sub)
    pos = layout(G)
    wmax = max(d["w"] for _, _, d in G.edges(data=True))
    for u, v, d in G.edges(data=True):
        span = abs(LAYER[u[1]] - LAYER[v[1]])
        nx.draw_networkx_edges(G, pos, edgelist=[(u, v)], ax=ax, width=0.8 + 5 * d["w"] / wmax,
                               edge_color=EDGE_COLOR[d["rel"]], alpha=0.75, arrows=True, arrowsize=10,
                               connectionstyle=f"arc3,rad={rad_for(span, pos[u][1], pos[v][1]):.3f}", node_size=1500)
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=[NODE_COLOR[n[1]] for n in G.nodes], node_size=1500,
                           edgecolors="black", linewidths=0.6)
    nx.draw_networkx_labels(G, pos, ax=ax, labels={n: textwrap.fill(n[0], 16) for n in G.nodes}, font_size=6.5)
    ax.set_title(title, fontsize=9.5)
    ax.axis("off")


def legend_handles():
    h = [Line2D([0], [0], marker="o", color="w", markerfacecolor=c, markeredgecolor="black", ms=9, label=l)
         for l, c in (("compound", NODE_COLOR["compound"]), ("bound gene", NODE_COLOR["bgene"]),
                      ("pathway", NODE_COLOR["pathway"]), ("disease gene", NODE_COLOR["dgene"]),
                      ("disease", NODE_COLOR["disease"]))]
    h += [Line2D([0], [0], color=c, lw=2, label=EDGE_LABEL[k]) for k, c in EDGE_COLOR.items()]
    return h


def path_title(r: pd.Series, per_mp: int, hidden: bool = False) -> str:
    note = "; disease-named pathways hidden" if hidden else ""
    return (f"{r['compound']}  [{r['tier']}; DWPC rank {int(r['dwpc_rank'])}]\n"
            f"binds {int(r['n_bound_genes'])} genes; top {per_mp} paths per metapath shown of "
            f"{int(r['n_paths_CbGaD'])} / {int(r['n_paths_CbGiGaD'])} / {int(r['n_paths_CbGpPWpGaD'])} "
            f"(CbGaD / CbGiGaD / CbGpPWpGaD){note}")


def drop_disease_pathway_paths(paths: pd.DataFrame, hide: set) -> pd.DataFrame:
    if not hide:
        return paths
    bad = paths["metapath"].eq("CbGpPWpGaD") & paths["path"].apply(
        lambda s: any(f" -in-> {h} <-in- " in s for h in hide))
    return paths[~bad]


def fig_paths(slug: str, summ: pd.DataFrame, paths: pd.DataFrame, hide: set, per_mp: int = 3) -> None:
    panels = []
    for _, r in summ.iterrows():
        cp = paths[paths["compound"] == r["compound"]]
        vis = drop_disease_pathway_paths(cp, hide)
        sub = vis.sort_values("share_in_metapath", ascending=False).groupby("metapath").head(per_mp)
        if sub.empty:
            continue
        title = path_title(r, per_mp, hidden=len(vis) < len(cp))
        fig, ax = plt.subplots(figsize=(11.5, 7))
        draw_candidate(ax, sub, title)
        ax.legend(handles=legend_handles(), loc="lower center", ncol=9, fontsize=7, frameon=False,
                  bbox_to_anchor=(0.5, -0.05))
        fig.tight_layout()
        safe = re.sub(r"[^A-Za-z0-9]+", "_", r["compound"])
        fig.savefig(OUT_FIG / f"path_{slug}_{safe}.png", dpi=200)
        plt.close(fig)
        panels.append((title, sub))
    if panels:
        top = panels[:4]
        fig, axes = plt.subplots(2, 2, figsize=(22, 15))
        for ax, (title, sub) in zip(axes.ravel(), top):
            draw_candidate(ax, sub, title)
        for ax in axes.ravel()[len(top):]:
            ax.axis("off")
        fig.legend(handles=legend_handles(), loc="lower center", ncol=9, fontsize=9, frameon=False)
        fig.suptitle(f"Supporting paths for top shortlisted candidates - {slug.capitalize()}'s "
                     "(edge width ~ share of the metapath's weight)", fontsize=13)
        fig.tight_layout(rect=(0, 0.03, 1, 0.97))
        fig.savefig(OUT_FIG / f"paths_{slug}_panel.png", dpi=170)
        plt.close(fig)
        print(f"  path_{slug}_*.png ({len(panels)}), paths_{slug}_panel.png")


def main() -> None:
    print("Writing figures to", OUT_FIG)
    nodes0, edges0 = load_processed()
    hide = set(build_matrices(nodes0, edges0, drop_disease_pathways=True).dropped_pathways or [])
    fig_benchmark()
    fig_damping()
    fig_hub_audit()
    all_rows = []
    for slug in DISEASES:
        summ = read(f"explain_summary_{slug}.csv")
        paths = read(f"paths_{slug}.csv")
        cand = read(f"candidates_{slug}.csv")
        if summ is None or paths is None or cand is None:
            continue
        summ = tier(summ)
        summ["disease"] = slug
        all_rows.append(summ)
        fig_evidence(slug, summ, cand)
        fig_paths(slug, summ, paths, hide)
    if all_rows:
        fig_ablation(pd.concat(all_rows, ignore_index=True))


if __name__ == "__main__":
    main()
