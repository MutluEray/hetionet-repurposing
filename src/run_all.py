"""Run the full pipeline end to end (skips steps whose outputs you do not need with --from / --to).

Order: 01 download -> 02 load -> 03 schema -> 04 RDF -> 05 SPARQL demo -> 06 scoring -> 07 evaluation ->
07b pathway check -> 08 permutation null -> 09 consensus -> 10 explanations -> 11 figures ->
12 interactive demo.
(03b rdflib timing is optional and not part of the default run.)

Usage:
    python3 src/run_all.py                 # everything (about 10 minutes; the null dominates)
    python3 src/run_all.py --n-perm 20     # quicker null
    python3 src/run_all.py --from 08       # resume from a step
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

SRC = Path(__file__).resolve().parent
STEPS = [
    ("01", "01_download.py", []),
    ("02", "02_load_hetionet.py", []),
    ("03", "03_schema_report.py", []),
    ("04", "04_build_rdf.py", []),
    ("05", "05_sparql_queries.py", []),
    ("06", "06_metapath_scoring.py", []),
    ("07", "07_evaluate.py", []),
    ("07b", "07b_pathway_leakage_check.py", []),
    ("08", "08_permutation_null.py", ["--n-perm", "{n_perm}"]),
    ("09", "09_consensus_candidates.py", []),
    ("10", "10_explain_paths.py", []),
    ("11", "11_figures.py", []),
    ("12", "12_interactive_demo.py", []),
]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--from", dest="start", default="01")
    ap.add_argument("--to", dest="stop", default="12")
    ap.add_argument("--n-perm", default="100")
    args = ap.parse_args()
    keys = [k for k, _, _ in STEPS]
    i0, i1 = keys.index(args.start), keys.index(args.stop)
    for key, script, extra in STEPS[i0:i1 + 1]:
        cmd = [sys.executable, str(SRC / script)] + [a.format(n_perm=args.n_perm) for a in extra]
        print(f"\n### {key}: {' '.join(cmd[1:])}", flush=True)
        t = time.time()
        if subprocess.run(cmd).returncode != 0:
            sys.exit(f"Step {key} failed.")
        print(f"### {key} done in {time.time() - t:.0f}s", flush=True)


if __name__ == "__main__":
    main()
