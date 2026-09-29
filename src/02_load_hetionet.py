"""Day 1: load Hetionet TSV/SIF into pandas, integer-index it, save as parquet.

Output (data/processed/):
    nodes.parquet : id, name, kind, ident, idx
    edges.parquet : metaedge, src_idx (int32), dst_idx (int32)
Edge endpoints are stored as integer indices (memory-light; ready for the
sparse-matrix scoring on Day 3). Use nodes.id[idx] to recover string IDs.
"""
from __future__ import annotations

import time

import numpy as np
import pandas as pd

from utils import (EDGES_PQ, EDGES_SIF, EXPECTED_EDGES, EXPECTED_NODES,
                   NODES_PQ, NODES_TSV)


def main() -> None:
    t0 = time.time()

    # ---- nodes ----------------------------------------------------------
    nodes = pd.read_csv(NODES_TSV, sep="\t", dtype=str)
    assert nodes.shape[1] == 3, f"Unexpected node columns: {list(nodes.columns)}"
    nodes.columns = ["id", "name", "kind"]
    nodes["ident"] = nodes["id"].str.split("::", n=1).str[1]
    prefix_mismatch = (nodes["id"].str.split("::", n=1).str[0] != nodes["kind"]).sum()
    nodes["kind"] = nodes["kind"].astype("category")
    nodes["idx"] = np.arange(len(nodes), dtype=np.int32)
    print(f"Nodes: {len(nodes):,} (expected {EXPECTED_NODES:,}); "
          f"id-prefix/kind mismatches: {prefix_mismatch}")

    # ---- edges ----------------------------------------------------------
    edges_raw = pd.read_csv(EDGES_SIF, sep="\t", dtype=str, compression="gzip")
    assert edges_raw.shape[1] == 3, f"Unexpected edge columns: {list(edges_raw.columns)}"
    edges_raw.columns = ["source", "metaedge", "target"]
    print(f"Edges: {len(edges_raw):,} (expected {EXPECTED_EDGES:,})")

    id_to_idx = pd.Series(nodes["idx"].to_numpy(), index=nodes["id"])
    src = edges_raw["source"].map(id_to_idx)
    dst = edges_raw["target"].map(id_to_idx)
    missing = int(src.isna().sum() + dst.isna().sum())
    print(f"Edge endpoints not found in node table: {missing}")
    if missing:
        raise SystemExit("Endpoint mismatch - inspect the data before continuing.")

    edges = pd.DataFrame({
        "metaedge": edges_raw["metaedge"].astype("category"),
        "src_idx": src.astype(np.int32).to_numpy(),
        "dst_idx": dst.astype(np.int32).to_numpy(),
    })
    del edges_raw

    dup = int(edges.duplicated(["src_idx", "metaedge", "dst_idx"]).sum())
    loops = int((edges["src_idx"] == edges["dst_idx"]).sum())
    print(f"Exact duplicate edges: {dup:,}; self-loops: {loops:,}")

    # ---- summaries ------------------------------------------------------
    print("\nNode counts by kind:")
    print(nodes["kind"].value_counts().to_string())

    print("\nEdge counts by metaedge (directed ones contain '>' or '<'):")
    print(edges["metaedge"].value_counts().to_string())

    # ---- save -----------------------------------------------------------
    nodes.to_parquet(NODES_PQ, index=False)
    edges.to_parquet(EDGES_PQ, index=False)

    mem_n = nodes.memory_usage(deep=True).sum() / 1e6
    mem_e = edges.memory_usage(deep=True).sum() / 1e6
    print(f"\nIn-memory size: nodes {mem_n:,.1f} MB, edges {mem_e:,.1f} MB")
    print(f"Saved {NODES_PQ.name}, {EDGES_PQ.name} to {NODES_PQ.parent}  ({time.time() - t0:.1f}s)")


if __name__ == "__main__":
    main()
