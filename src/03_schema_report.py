"""Day 1: schema report, degree analysis, and AD/PD neighborhood summary.

Writes to outputs/tables/ and outputs/figures/:
    node_counts.csv, metaedges.csv, degree_by_kind.csv, top_hubs_by_kind.csv,
    diseases_all.csv, <slug>_neighborhood.csv, <slug>_known_drugs.csv,
    scope_estimate.csv
    metagraph.png, degree_distributions.png
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch

from utils import OUT_FIG, OUT_TAB, find_nodes, load_processed, node_degree

TARGETS = [("alzheimer", "alzheimer"), ("parkinson", "parkinson")]  # (name query, file slug)


# --------------------------------------------------------------------------
def metaedge_table(nodes, edges):
    codes = nodes["kind"].cat.codes.to_numpy()
    cats = nodes["kind"].cat.categories
    first = edges.drop_duplicates("metaedge")[["metaedge", "src_idx", "dst_idx"]]
    kinds = pd.DataFrame({
        "metaedge": first["metaedge"].astype(str).to_numpy(),
        "src_kind": cats[codes[first["src_idx"].to_numpy()]],
        "dst_kind": cats[codes[first["dst_idx"].to_numpy()]],
    })
    g = edges.groupby("metaedge", observed=True)
    stats = pd.DataFrame({
        "n_edges": g.size(),
        "n_unique_src": g["src_idx"].nunique(),
        "n_unique_dst": g["dst_idx"].nunique(),
    }).reset_index()
    stats["metaedge"] = stats["metaedge"].astype(str)
    out = kinds.merge(stats, on="metaedge").sort_values("n_edges", ascending=False)
    out["directed"] = out["metaedge"].str.contains("[<>]")
    return out.reset_index(drop=True)


def plot_metagraph(nodes, met, path):
    kinds = list(nodes["kind"].cat.categories)
    counts = nodes["kind"].value_counts()
    n = len(kinds)
    ang = np.linspace(0, 2 * np.pi, n, endpoint=False) + np.pi / 2
    pos = {k: (np.cos(a), np.sin(a)) for k, a in zip(kinds, ang)}

    pair = {}
    for _, r in met.iterrows():
        key = tuple(sorted((r.src_kind, r.dst_kind)))
        pair[key] = pair.get(key, 0) + r.n_edges

    fig, ax = plt.subplots(figsize=(10, 10))
    for (a, b), c in pair.items():
        lw = 0.4 + 0.9 * np.log10(max(c, 10))
        if a == b:
            x, y = pos[a]
            ax.add_patch(Circle((x * 1.12, y * 1.12), 0.07, fill=False, lw=lw, alpha=0.5, color="tab:blue"))
        else:
            ax.add_patch(FancyArrowPatch(pos[a], pos[b], arrowstyle="-", lw=lw, alpha=0.4,
                                         color="tab:blue", connectionstyle="arc3,rad=0.08"))
    for k in kinds:
        x, y = pos[k]
        ax.scatter([x], [y], s=900, color="white", edgecolor="black", zorder=3)
        ax.text(x * 1.0, y * 1.0, f"{k}\n{counts[k]:,}", ha="center", va="center", fontsize=8, zorder=4)
    ax.set_xlim(-1.4, 1.4)
    ax.set_ylim(-1.4, 1.4)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Hetionet v1.0 metagraph (node type counts; line width ~ log edge count)")
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_degrees(nodes, path):
    kinds = list(nodes["kind"].cat.categories)
    fig, axes = plt.subplots(3, 4, figsize=(16, 10))
    for ax, k in zip(axes.ravel(), kinds):
        d = nodes.loc[nodes["kind"] == k, "degree"].to_numpy()
        d = d[d > 0]
        bins = np.logspace(0, np.log10(max(d.max(), 2)), 30)
        ax.hist(d, bins=bins)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(f"{k} (n={len(d):,})", fontsize=9)
        ax.set_xlabel("degree")
    for ax in axes.ravel()[len(kinds):]:
        ax.axis("off")
    fig.suptitle("Degree distributions by node type (log-log)")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def neighbor_table(nodes, edges, d_idx):
    src = edges["src_idx"].to_numpy()
    dst = edges["dst_idx"].to_numpy()
    out = edges[src == d_idx][["metaedge", "dst_idx"]].rename(columns={"dst_idx": "nbr_idx"})
    out["direction"] = "out"
    inn = edges[dst == d_idx][["metaedge", "src_idx"]].rename(columns={"src_idx": "nbr_idx"})
    inn["direction"] = "in"
    t = pd.concat([out, inn], ignore_index=True)
    t["metaedge"] = t["metaedge"].astype(str)
    t["nbr_id"] = nodes["id"].to_numpy()[t["nbr_idx"]]
    t["nbr_name"] = nodes["name"].to_numpy()[t["nbr_idx"]]
    t["nbr_kind"] = nodes["kind"].astype(str).to_numpy()[t["nbr_idx"]]
    return t


def main() -> None:
    nodes, edges = load_processed()
    nodes["degree"] = node_degree(nodes, edges)

    # ---- global tables --------------------------------------------------
    nodes["kind"].value_counts().rename_axis("kind").reset_index(name="n_nodes") \
        .to_csv(OUT_TAB / "node_counts.csv", index=False)

    met = metaedge_table(nodes, edges)
    met.to_csv(OUT_TAB / "metaedges.csv", index=False)
    print("Metaedges:")
    print(met.to_string(index=False))

    deg_stats = nodes.groupby("kind", observed=True)["degree"].agg(
        ["count", "mean", "median", "max"]).round(1)
    deg_stats.to_csv(OUT_TAB / "degree_by_kind.csv")
    print("\nDegree by node kind:")
    print(deg_stats.to_string())

    hubs = (nodes.sort_values("degree", ascending=False)
            .groupby("kind", observed=True).head(10)[["kind", "id", "name", "degree"]]
            .sort_values(["kind", "degree"], ascending=[True, False]))
    hubs.to_csv(OUT_TAB / "top_hubs_by_kind.csv", index=False)

    plot_metagraph(nodes, met, OUT_FIG / "metagraph.png")
    plot_degrees(nodes, OUT_FIG / "degree_distributions.png")

    # ---- all diseases, with treatment/association counts ------------------
    is_d = (nodes["kind"].astype(str) == "Disease").to_numpy()
    s, d = edges["src_idx"].to_numpy(), edges["dst_idx"].to_numpy()
    ms, md = is_d[s], is_d[d]
    touch = pd.concat([
        pd.DataFrame({"idx": s[ms], "metaedge": edges["metaedge"].astype(str).to_numpy()[ms]}),
        pd.DataFrame({"idx": d[md], "metaedge": edges["metaedge"].astype(str).to_numpy()[md]}),
    ])
    piv = touch.groupby(["idx", "metaedge"]).size().unstack(fill_value=0)
    dis = nodes[is_d][["idx", "id", "name", "degree"]].merge(piv, left_on="idx", right_index=True, how="left").fillna(0)
    dis.sort_values("degree", ascending=False).to_csv(OUT_TAB / "diseases_all.csv", index=False)
    print(f"\n{len(dis)} disease nodes saved to diseases_all.csv "
          f"(check CtD column for how many known treatments each has).")

    # ---- AD / PD neighborhoods ---------------------------------------------
    scope_rows = []
    for query, slug in TARGETS:
        hits = find_nodes(nodes, query, kind="Disease")
        print(f"\n=== {slug.upper()} ===")
        if hits.empty:
            print(f"No Disease node matching '{query}'. See diseases_all.csv.")
            continue
        print("Matches:")
        print(hits[["id", "name", "degree"]].to_string(index=False))
        row = hits.sort_values("degree", ascending=False).iloc[0]
        d_idx = int(row["idx"])
        print(f"Using: {row['id']}  ({row['name']})")

        nb = neighbor_table(nodes, edges, d_idx)
        summ = (nb.groupby(["metaedge", "direction", "nbr_kind"]).size()
                .reset_index(name="n").sort_values("n", ascending=False))
        summ.to_csv(OUT_TAB / f"{slug}_neighborhood.csv", index=False)
        print(summ.to_string(index=False))

        drugs = nb[nb["metaedge"].isin(["CtD", "CpD"])][["metaedge", "nbr_id", "nbr_name"]] \
            .rename(columns={"nbr_id": "compound_id", "nbr_name": "compound"})
        drugs.to_csv(OUT_TAB / f"{slug}_known_drugs.csv", index=False)
        print(f"\nKnown drugs in graph (CtD = treats, CpD = palliates): {len(drugs)}")
        print(drugs.to_string(index=False))

        # scope estimate for the Day-2 subgraph
        genes = set(nb[nb["metaedge"].isin(["DaG", "DuG", "DdG"])]["nbr_idx"])
        g_arr = np.fromiter(genes, dtype=np.int64) if genes else np.array([], dtype=np.int64)
        me = edges["metaedge"].astype(str).to_numpy()
        cbg = (me == "CbG") & np.isin(d, g_arr)
        gpw = (me == "GpPW") & np.isin(s, g_arr)
        gig = (me == "GiG") & (np.isin(s, g_arr) | np.isin(d, g_arr))
        scope_rows.append({
            "disease": row["name"],
            "genes_linked_DaG_DuG_DdG": len(genes),
            "CbG_edges_into_those_genes": int(cbg.sum()),
            "compounds_binding_those_genes": int(len(np.unique(s[cbg]))),
            "GpPW_edges_from_those_genes": int(gpw.sum()),
            "pathways_touched": int(len(np.unique(d[gpw]))),
            "GiG_edges_touching_those_genes": int(gig.sum()),
        })

    if scope_rows:
        sc = pd.DataFrame(scope_rows)
        sc.to_csv(OUT_TAB / "scope_estimate.csv", index=False)
        print("\nSubgraph scope estimate (for Day 2 scoping decisions):")
        print(sc.T.to_string())


if __name__ == "__main__":
    main()
