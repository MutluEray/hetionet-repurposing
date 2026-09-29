"""Day 1: extract the Hetionet files we need from the Zenodo bundle.

Expects the bundle (hetionet-v1.0.0.zip from https://zenodo.org/records/268568)
somewhere in data/raw/. Extracts only the TSV / JSON / describe files (skips the
big Neo4j dump), copies them to data/raw/hetionet/, and reports actual file sizes.

Usage:
    python src/01_download.py            # normal run
    python src/01_download.py --force    # re-extract even if files already exist
    python src/01_download.py --keep     # keep the temporary extraction folder
"""
from __future__ import annotations

import argparse
import shutil
import sys
import zipfile
from pathlib import Path

from utils import DATA_RAW, EDGES_SIF, HET_DIR, NODES_TSV

TMP = DATA_RAW / "_extracted"
WANT_DIRS = ("hetnet/tsv/", "hetnet/json/", "describe/")


def human(n: float) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:,.1f} {unit}"
        n /= 1024
    return f"{n:.1f} GB"


def find_zip() -> Path:
    zips = sorted(DATA_RAW.glob("*.zip"))
    if not zips:
        sys.exit(f"No .zip found in {DATA_RAW}. Put hetionet-v1.0.0.zip there and re-run.")
    preferred = [z for z in zips if "hetionet" in z.name.lower()]
    return (preferred or zips)[0]


def extract_needed(zip_path: Path) -> None:
    if TMP.exists():
        shutil.rmtree(TMP)
    TMP.mkdir(parents=True)
    with zipfile.ZipFile(zip_path) as zf:
        names = zf.namelist()
        print(f"Zip contains {len(names):,} entries; top level: "
              f"{sorted({n.split('/')[0] for n in names})[:5]}")
        wanted = [n for n in names
                  if not n.endswith("/") and any(d in n.lower() for d in WANT_DIRS)]
        if not wanted:
            print("No entries matched expected folders; extracting everything.")
            wanted = [n for n in names if not n.endswith("/")]
        for n in wanted:
            zf.extract(n, TMP)
        print(f"Extracted {len(wanted)} files to {TMP}")


def first(pattern: str) -> Path | None:
    hits = sorted(TMP.rglob(pattern))
    return hits[0] if hits else None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--keep", action="store_true")
    args = ap.parse_args()

    if NODES_TSV.exists() and EDGES_SIF.exists() and not args.force:
        print("Hetionet files already in place (use --force to redo).")
    else:
        zip_path = find_zip()
        print(f"Using {zip_path} ({human(zip_path.stat().st_size)})")
        extract_needed(zip_path)

        nodes = first("hetionet-v1.0-nodes.tsv")
        edges = first("hetionet-v1.0-edges.sif.gz")
        if nodes is None or edges is None:
            print("Could not find nodes.tsv / edges.sif.gz. Files seen containing 'hetionet' or 'tsv':")
            for p in sorted(TMP.rglob("*")):
                if p.is_file() and ("hetionet" in p.name.lower() or p.suffix == ".tsv"):
                    print("  ", p.relative_to(TMP))
            sys.exit(1)

        shutil.copy2(nodes, NODES_TSV)
        shutil.copy2(edges, EDGES_SIF)

        # Optional extras
        js = first("hetionet-v1.0.json.bz2")
        if js:
            shutil.copy2(js, HET_DIR / js.name)
        for name in ("metaedges.tsv", "metanodes.tsv"):
            p = first(name)
            if p:
                shutil.copy2(p, HET_DIR / name)

        if not args.keep:
            shutil.rmtree(TMP)

    print("\nFiles in data/raw/hetionet/ (actual sizes):")
    for p in sorted(HET_DIR.iterdir()):
        if p.is_file():
            print(f"  {p.name:<40} {human(p.stat().st_size):>12}")


if __name__ == "__main__":
    main()
