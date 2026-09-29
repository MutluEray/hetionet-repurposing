"""Day 2: SPARQL exploration of the Hetionet RDF graph with rdflib.

Loads outputs/rdf/hetionet_core.nt (build it with 04_build_rdf.py), checks the
RDF round-trip against the pandas edge counts, then runs anchored multi-hop
queries for one target disease:

    Q1  known drugs (CtD / CpD)
    Q2  CbGaD      compounds binding disease-associated genes
    Q3  CbGpPWpGaD compounds sharing a pathway with disease-associated genes  (heavier)
    Q4  CtDrD      compounds treating diseases that resemble the target (property path for symmetric DrD)
    Q5  signature reversal (CdGuD + CuGdD)
    Q6  explicit path retrieval for the top candidates (the explanation step)

Known drugs are excluded from Q2-Q5 (candidate view). Results go to outputs/tables/.

Usage:
    python3 src/05_sparql_queries.py
    python3 src/05_sparql_queries.py --disease parkinson
    python3 src/05_sparql_queries.py --skip-heavy          # skip Q3
    python3 src/05_sparql_queries.py --graph outputs/rdf/hetionet_all.nt

NOTE: keep queries anchored on a bound disease/compound. Unanchored multi-hop joins
over hubs explode (a 2-hop join over a 100k-edge sample already returned 340k rows).
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path
from urllib.parse import unquote

import pandas as pd
from rdflib import Graph

from utils import OUT_RDF, OUT_TAB, find_nodes, load_processed, node_uri

PFX = """\
PREFIX rel:  <http://het.io/rel/>
PREFIX rdf:  <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
"""

Q_TYPES = """
SELECT ?type (COUNT(?n) AS ?count) WHERE { ?n rdf:type ?type } GROUP BY ?type ORDER BY DESC(?count)
"""

Q_REL_COUNTS = """
SELECT ?p (COUNT(*) AS ?n) WHERE {
  ?s ?p ?o . FILTER(STRSTARTS(STR(?p), "http://het.io/rel/"))
} GROUP BY ?p
"""

Q1_KNOWN = """
SELECT ?relation ?compound WHERE {
  VALUES ?relation { rel:CtD rel:CpD }
  ?c ?relation %D% .
  ?c rdfs:label ?compound .
} ORDER BY ?relation ?compound
"""

Q2_CBGAD = """
SELECT ?c ?compound (COUNT(DISTINCT ?g) AS ?n_shared_genes)
       (GROUP_CONCAT(DISTINCT ?gname; separator=", ") AS ?genes)
WHERE {
  %D% rel:DaG ?g .
  ?c rel:CbG ?g .
  FILTER NOT EXISTS { ?c rel:CtD %D% }
  FILTER NOT EXISTS { ?c rel:CpD %D% }
  ?c rdfs:label ?compound .
  ?g rdfs:label ?gname .
}
GROUP BY ?c ?compound
ORDER BY DESC(?n_shared_genes)
LIMIT 25
"""

Q3_PATHWAY = """
SELECT ?c ?compound (COUNT(DISTINCT ?p) AS ?n_pathways) (COUNT(*) AS ?n_paths)
WHERE {
  %D% rel:DaG ?g2 .
  ?g2 rel:GpPW ?p .
  ?g1 rel:GpPW ?p .
  ?c rel:CbG ?g1 .
  FILTER NOT EXISTS { ?c rel:CtD %D% }
  FILTER NOT EXISTS { ?c rel:CpD %D% }
  ?c rdfs:label ?compound .
}
GROUP BY ?c ?compound
ORDER BY DESC(?n_pathways)
LIMIT 25
"""

Q4_SIMILAR = """
SELECT ?c ?compound ?similar_disease WHERE {
  %D% (rel:DrD|^rel:DrD) ?d2 .
  ?c rel:CtD ?d2 .
  FILTER NOT EXISTS { ?c rel:CtD %D% }
  FILTER NOT EXISTS { ?c rel:CpD %D% }
  ?c rdfs:label ?compound .
  ?d2 rdfs:label ?similar_disease .
} ORDER BY ?similar_disease ?compound
LIMIT 50
"""

Q5_REVERSAL = """
SELECT ?c ?compound (COUNT(DISTINCT ?g) AS ?n_reversed_genes) WHERE {
  { %D% rel:DuG ?g . ?c rel:CdG ?g }
  UNION
  { %D% rel:DdG ?g . ?c rel:CuG ?g }
  FILTER NOT EXISTS { ?c rel:CtD %D% }
  FILTER NOT EXISTS { ?c rel:CpD %D% }
  ?c rdfs:label ?compound .
}
GROUP BY ?c ?compound
ORDER BY DESC(?n_reversed_genes)
LIMIT 25
"""

Q6_PATH_GENES = """
SELECT ?compound ?gene WHERE {
  %C% rel:CbG ?g .
  %D% rel:DaG ?g .
  %C% rdfs:label ?compound .
  ?g rdfs:label ?gene .
} ORDER BY ?gene
"""

Q6_PATH_PATHWAY = """
SELECT ?compound ?bound_gene ?pathway ?disease_gene WHERE {
  %C% rel:CbG ?g1 .
  ?g1 rel:GpPW ?p .
  ?g2 rel:GpPW ?p .
  %D% rel:DaG ?g2 .
  %C% rdfs:label ?compound .
  ?g1 rdfs:label ?bound_gene .
  ?p rdfs:label ?pathway .
  ?g2 rdfs:label ?disease_gene .
  FILTER(?g1 != ?g2)
} LIMIT 20
"""


def run(g: Graph, title: str, query: str, csv_name: str | None = None, show: int = 15) -> pd.DataFrame:
    t = time.time()
    res = g.query(PFX + query)
    cols = [str(v) for v in res.vars]
    rows = [[("" if v is None else str(v)) for v in r] for r in res]
    df = pd.DataFrame(rows, columns=cols)
    for c in df.columns:
        conv = pd.to_numeric(df[c], errors="coerce")
        if len(df) and conv.notna().all():
            df[c] = conv
    print(f"\n--- {title}  ({time.time() - t:.2f}s, {len(df)} rows)")
    if len(df):
        disp = df.drop(columns=[c for c in ("c",) if c in df.columns]).head(show).copy()
        for c in disp.columns:
            if disp[c].dtype == object:
                disp[c] = disp[c].astype(str).str.slice(0, 70)
        print(disp.to_string(index=False))
    if csv_name:
        df.to_csv(OUT_TAB / csv_name, index=False)
    return df


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--graph", default=str(OUT_RDF / "hetionet_core.nt"))
    ap.add_argument("--disease", default="alzheimer", help="name substring of a Disease node")
    ap.add_argument("--skip-heavy", action="store_true")
    args = ap.parse_args()

    nodes, edges = load_processed()
    hits = find_nodes(nodes, args.disease, kind="Disease")
    if hits.empty:
        raise SystemExit(f"No Disease node matching '{args.disease}'")
    row = hits.iloc[0]
    slug = args.disease.lower().replace(" ", "_")
    D = f"<{node_uri(row['id'])}>"
    print(f"Target disease: {row['name']} ({row['id']})")

    gpath = Path(args.graph)
    if not gpath.exists():
        raise SystemExit(f"{gpath} not found. Run 04_build_rdf.py first.")
    t = time.time()
    g = Graph()
    g.parse(str(gpath), format="nt")
    print(f"Loaded {len(g):,} triples in {time.time() - t:.1f}s")

    # ---- census + round-trip check -------------------------------------------
    run(g, "Node counts by type", Q_TYPES)
    res = g.query(PFX + Q_REL_COUNTS)
    rdf_counts = {unquote(str(r[0]).replace("http://het.io/rel/", "")): int(r[1]) for r in res}
    pd_counts = edges["metaedge"].astype(str).value_counts().to_dict()
    bad = {m: (n, pd_counts.get(m)) for m, n in rdf_counts.items() if pd_counts.get(m) != n}
    print("\nRound-trip check (RDF vs pandas edge counts per relation): "
          + ("OK" if not bad else f"MISMATCH {bad}"))

    # ---- disease queries ------------------------------------------------------
    sub = lambda q: q.replace("%D%", D)
    run(g, "Q1 known drugs (CtD/CpD)", sub(Q1_KNOWN), f"sparql_{slug}_q1_known.csv", show=40)
    q2 = run(g, "Q2 CbGaD: compounds binding disease-associated genes", sub(Q2_CBGAD), f"sparql_{slug}_q2_cbgad.csv")
    q3 = None
    if not args.skip_heavy:
        q3 = run(g, "Q3 CbGpPWpGaD: shared-pathway compounds", sub(Q3_PATHWAY), f"sparql_{slug}_q3_pathway.csv")
    run(g, "Q4 CtDrD: compounds treating similar diseases", sub(Q4_SIMILAR), f"sparql_{slug}_q4_similar.csv", show=30)
    run(g, "Q5 signature reversal (CdGuD + CuGdD)", sub(Q5_REVERSAL), f"sparql_{slug}_q5_reversal.csv")

    # ---- explanation step: explicit paths for top candidates -------------------
    if len(q2):
        C = f"<{q2.iloc[0]['c']}>"
        run(g, f"Q6a genes on the CbGaD path of top Q2 candidate ({q2.iloc[0]['compound']})",
            sub(Q6_PATH_GENES).replace("%C%", C), f"sparql_{slug}_q6_paths_cbgad.csv", show=30)
    if q3 is not None and len(q3):
        C = f"<{q3.iloc[0]['c']}>"
        run(g, f"Q6b pathway paths of top Q3 candidate ({q3.iloc[0]['compound']})",
            sub(Q6_PATH_PATHWAY).replace("%C%", C), f"sparql_{slug}_q6_paths_pathway.csv", show=20)


if __name__ == "__main__":
    main()
