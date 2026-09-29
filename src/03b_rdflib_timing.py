"""Day 1: empirical rdflib go/no-go test.

Builds an rdflib graph from a random sample of Hetionet edges, times build /
serialize / parse / two SPARQL queries, measures memory, and extrapolates
(crudely, linearly) to the full 2.25M-edge graph.

Usage:
    python src/03b_rdflib_timing.py            # 100k-edge sample
    python src/03b_rdflib_timing.py --n 300000
"""
from __future__ import annotations

import argparse
import gc
import os
import subprocess
import time

from rdflib import Graph, URIRef

from utils import DATA_PROC, edge_uri, load_processed, node_uri


def rss_mb() -> float:
    """Current resident memory of this process in MB (via `ps`; no extra dependency)."""
    out = subprocess.check_output(["ps", "-o", "rss=", "-p", str(os.getpid())])
    return int(out.strip()) / 1024.0


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=100_000)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    nodes, edges = load_processed()
    total = len(edges)
    sample = edges.sample(n=min(args.n, total), random_state=args.seed)
    ids = nodes["id"].to_numpy()
    s_ids = ids[sample["src_idx"].to_numpy()]
    o_ids = ids[sample["dst_idx"].to_numpy()]
    mets = sample["metaedge"].astype(str).to_numpy()
    n = len(sample)
    scale = total / n
    print(f"Sample: {n:,} of {total:,} edges (scale factor to full graph: {scale:.1f}x)")

    cache: dict[str, URIRef] = {}

    def U(node_id: str) -> URIRef:
        u = cache.get(node_id)
        if u is None:
            u = cache[node_id] = URIRef(node_uri(node_id))
        return u

    pred = {m: URIRef(edge_uri(m)) for m in set(mets)}

    # ---- build ----------------------------------------------------------
    gc.collect()
    rss0 = rss_mb()
    t = time.time()
    g = Graph()
    for s, m, o in zip(s_ids, mets, o_ids):
        g.add((U(s), pred[m], U(o)))
    build_t = time.time() - t
    build_mem = rss_mb() - rss0
    print(f"\nBuild:     {build_t:6.1f}s   +{build_mem:,.0f} MB   ({len(g):,} triples)")

    # ---- serialize ------------------------------------------------------
    nt = DATA_PROC / "rdflib_timing_sample.nt"
    t = time.time()
    g.serialize(destination=str(nt), format="nt")
    ser_t = time.time() - t
    print(f"Serialize: {ser_t:6.1f}s   file {nt.stat().st_size / 1e6:,.1f} MB")

    # ---- parse ----------------------------------------------------------
    del g
    cache.clear()
    gc.collect()
    rss1 = rss_mb()
    t = time.time()
    g2 = Graph()
    g2.parse(str(nt), format="nt")
    parse_t = time.time() - t
    parse_mem = rss_mb() - rss1
    print(f"Parse:     {parse_t:6.1f}s   +{parse_mem:,.0f} MB")

    # ---- queries --------------------------------------------------------
    q_count = "SELECT (COUNT(*) AS ?n) WHERE { ?s ?p ?o }"
    q_join = "SELECT (COUNT(*) AS ?n) WHERE { ?a ?p1 ?b . ?b ?p2 ?c }"
    for label, q in (("COUNT all triples", q_count), ("2-hop join count", q_join)):
        t = time.time()
        res = list(g2.query(q))
        print(f"Query [{label}]: {time.time() - t:6.2f}s  -> {res[0][0]}")

    # ---- extrapolation ----------------------------------------------------
    print("\nCrude linear extrapolation to the full graph (joins scale worse than linearly):")
    print(f"  build : ~{build_t * scale / 60:,.1f} min, ~{build_mem * scale / 1000:,.1f} GB")
    print(f"  parse : ~{parse_t * scale / 60:,.1f} min, ~{parse_mem * scale / 1000:,.1f} GB")
    print("Decision rule: if these are uncomfortable for your Mac, keep rdflib on the scoped subgraph only.")


if __name__ == "__main__":
    main()
