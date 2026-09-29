"""Shared helpers for the Hetionet drug-repurposing pipeline.

Paths, processed-data loaders, node lookup, degree computation and the
URI scheme used when exporting to RDF (reused on Day 2).
"""
from __future__ import annotations

from pathlib import Path
from urllib.parse import quote

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parent.parent
DATA_RAW = ROOT / "data" / "raw"
DATA_PROC = ROOT / "data" / "processed"
HET_DIR = DATA_RAW / "hetionet"          # clean copies of the files we use
OUT_FIG = ROOT / "outputs" / "figures"
OUT_TAB = ROOT / "outputs" / "tables"
OUT_RDF = ROOT / "outputs" / "rdf"

for _d in (DATA_RAW, DATA_PROC, HET_DIR, OUT_FIG, OUT_TAB, OUT_RDF):
    _d.mkdir(parents=True, exist_ok=True)

NODES_TSV = HET_DIR / "hetionet-v1.0-nodes.tsv"
EDGES_SIF = HET_DIR / "hetionet-v1.0-edges.sif.gz"
NODES_PQ = DATA_PROC / "nodes.parquet"
EDGES_PQ = DATA_PROC / "edges.parquet"

EXPECTED_NODES = 47_031
EXPECTED_EDGES = 2_250_197

# --------------------------------------------------------------------------
# RDF URI scheme (used from Day 2 on)
# --------------------------------------------------------------------------
HET_NODE_NS = "http://het.io/node/"
HET_EDGE_NS = "http://het.io/rel/"


HET_CLASS_NS = "http://het.io/class/"
HET_ONT_NS = "http://het.io/ns/"

# Human-readable labels for the 24 Hetionet metaedges.
METAEDGE_LABELS = {
    "AdG": "Anatomy downregulates Gene", "AeG": "Anatomy expresses Gene",
    "AuG": "Anatomy upregulates Gene", "CbG": "Compound binds Gene",
    "CcSE": "Compound causes Side Effect", "CdG": "Compound downregulates Gene",
    "CpD": "Compound palliates Disease", "CrC": "Compound resembles Compound",
    "CtD": "Compound treats Disease", "CuG": "Compound upregulates Gene",
    "DaG": "Disease associates Gene", "DdG": "Disease downregulates Gene",
    "DlA": "Disease localizes Anatomy", "DpS": "Disease presents Symptom",
    "DrD": "Disease resembles Disease", "DuG": "Disease upregulates Gene",
    "GcG": "Gene covaries Gene", "GiG": "Gene interacts Gene",
    "Gr>G": "Gene regulates Gene", "GpBP": "Gene participates Biological Process",
    "GpCC": "Gene participates Cellular Component", "GpMF": "Gene participates Molecular Function",
    "GpPW": "Gene participates Pathway", "PCiC": "Pharmacologic Class includes Compound",
}

# Mechanistic "core": relations used by the repurposing metapaths (drops the large
# annotation/hub relations: GpBP, AeG, AdG, AuG, GpMF, GpCC, CcSE, Gr>G, GcG, DlA, DpS).
CORE_RELATIONS = ["CbG", "CuG", "CdG", "CrC", "PCiC", "CtD", "CpD",
                  "DaG", "DuG", "DdG", "DrD", "GpPW", "GiG"]


def class_uri(kind: str) -> str:
    """'Biological Process' -> 'http://het.io/class/BiologicalProcess'."""
    return HET_CLASS_NS + quote(kind.replace(" ", ""), safe="")


def node_uri(node_id: str) -> str:
    """'Gene::3149' -> 'http://het.io/node/Gene/3149' (percent-encoded)."""
    kind, ident = node_id.split("::", 1)
    return f"{HET_NODE_NS}{quote(kind.replace(' ', ''), safe='')}/{quote(ident, safe='')}"


def edge_uri(metaedge: str) -> str:
    """'Gr>G' -> 'http://het.io/rel/Gr%3EG'."""
    return HET_EDGE_NS + quote(metaedge, safe="")


# --------------------------------------------------------------------------
# Loading
# --------------------------------------------------------------------------
def load_processed() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (nodes, edges) written by 02_load_hetionet.py.

    nodes: id, name, kind (category), ident, idx
    edges: metaedge (category), src_idx (int32), dst_idx (int32)
    """
    if not NODES_PQ.exists() or not EDGES_PQ.exists():
        raise FileNotFoundError(
            "Processed parquet files not found. Run 01_download.py and 02_load_hetionet.py first."
        )
    nodes = pd.read_parquet(NODES_PQ)
    edges = pd.read_parquet(EDGES_PQ)
    nodes["kind"] = nodes["kind"].astype("category")
    edges["metaedge"] = edges["metaedge"].astype("category")
    return nodes, edges


def node_degree(nodes: pd.DataFrame, edges: pd.DataFrame) -> np.ndarray:
    """Total degree (in + out) per node, aligned with nodes.idx."""
    n = len(nodes)
    return (
        np.bincount(edges["src_idx"].to_numpy(), minlength=n)
        + np.bincount(edges["dst_idx"].to_numpy(), minlength=n)
    )


def find_nodes(nodes: pd.DataFrame, query: str, kind: str | None = None) -> pd.DataFrame:
    """Case-insensitive substring search on node names, optionally restricted to a kind."""
    mask = nodes["name"].str.contains(query, case=False, regex=False, na=False)
    if kind is not None:
        mask &= nodes["kind"].astype(str) == kind
    return nodes[mask]
