"""Day 2: export Hetionet (or a relation-scoped subset) to RDF N-Triples.

Scoping is done by *relation type*, not by disease neighborhood: the default
"core" preset keeps the 13 mechanistic relations used by the repurposing
metapaths (~320k edges) and drops the large annotation/hub relations
(GpBP, AeG, AdG, AuG, GpMF, GpCC, CcSE, Gr>G, GcG, DlA, DpS). Metapath results on
the core graph are therefore identical to those on the full graph.

Each node gets rdf:type (its kind) and rdfs:label (its name); each relation gets
an rdfs:label plus source/target kind, so the graph is self-describing.

Usage:
    python3 src/04_build_rdf.py                          # core preset -> outputs/rdf/hetionet_core.nt
    python3 src/04_build_rdf.py --relations all          # full graph (2.25M edges, ~230 MB file)
    python3 src/04_build_rdf.py --relations CbG,DaG,GpPW --out mini
"""
from __future__ import annotations

import argparse
import time

import numpy as np
from rdflib import Literal

from utils import (CORE_RELATIONS, HET_ONT_NS, METAEDGE_LABELS, OUT_RDF,
                   class_uri, edge_uri, load_processed, node_uri)

RDF_TYPE = "<http://www.w3.org/1999/02/22-rdf-syntax-ns#type>"
RDFS_LABEL = "<http://www.w3.org/2000/01/rdf-schema#label>"
SRC_KIND = f"<{HET_ONT_NS}sourceKind>"
DST_KIND = f"<{HET_ONT_NS}targetKind>"


def resolve_relations(arg: str, all_rels: list[str]) -> list[str]:
    if arg == "core":
        rels = list(CORE_RELATIONS)
    elif arg == "all":
        rels = list(all_rels)
    else:
        rels = [r.strip() for r in arg.split(",") if r.strip()]
    unknown = [r for r in rels if r not in set(all_rels)]
    if unknown:
        raise SystemExit(f"Unknown relation(s): {unknown}. Available: {all_rels}")
    return rels


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--relations", default="core",
                    help="'core', 'all', or comma-separated metaedge abbreviations")
    ap.add_argument("--out", default=None, help="output file stem")
    args = ap.parse_args()

    t0 = time.time()
    nodes, edges = load_processed()
    all_rels = [str(x) for x in edges["metaedge"].cat.categories]
    rels = resolve_relations(args.relations, all_rels)
    preset = args.relations if args.relations in ("core", "all") else "custom"
    stem = args.out or f"hetionet_{preset}"
    path = OUT_RDF / f"{stem}.nt"

    met_all = edges["metaedge"].astype(str).to_numpy()
    keep = np.isin(met_all, rels)
    src = edges["src_idx"].to_numpy()[keep].tolist()
    dst = edges["dst_idx"].to_numpy()[keep].tolist()
    met = met_all[keep].tolist()
    print(f"Relations ({len(rels)}): {', '.join(rels)}")
    print(f"Selected {len(met):,} of {len(edges):,} edges")

    ids = nodes["id"].to_numpy()
    names = nodes["name"].to_numpy()
    kinds = nodes["kind"].astype(str).to_numpy()
    used = sorted(set(src) | set(dst))
    uri = {i: f"<{node_uri(ids[i])}>" for i in used}
    pred = {m: f"<{edge_uri(m)}>" for m in rels}

    n_triples = 0
    with open(path, "w", encoding="utf-8") as f:
        def w(s: str, p: str, o: str) -> None:
            nonlocal n_triples
            f.write(f"{s} {p} {o} .\n")
            n_triples += 1

        # relation metadata
        first_seen: dict[str, int] = {}
        for k, m in enumerate(met):
            if m not in first_seen:
                first_seen[m] = k
        for m in rels:
            w(pred[m], RDFS_LABEL, Literal(METAEDGE_LABELS.get(m, m)).n3())
            if m in first_seen:
                k = first_seen[m]
                w(pred[m], SRC_KIND, f"<{class_uri(kinds[src[k]])}>")
                w(pred[m], DST_KIND, f"<{class_uri(kinds[dst[k]])}>")

        # nodes: type + label
        for i in used:
            w(uri[i], RDF_TYPE, f"<{class_uri(kinds[i])}>")
            w(uri[i], RDFS_LABEL, Literal(str(names[i])).n3())

        # edges
        for s, m, o in zip(src, met, dst):
            w(uri[s], pred[m], uri[o])

    print(f"\nNodes written: {len(used):,}")
    print(f"Triples written: {n_triples:,}  (edges {len(met):,} + node type/label + relation metadata)")
    print(f"File: {path}  ({path.stat().st_size / 1e6:,.1f} MB)  in {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
