"""Day 7: temporal validation against ClinicalTrials.gov (implements docs/validation_plan.md).

Question: does the frozen consensus score rank compounds that LATER entered phase 2/3 trials for Alzheimer's or
Parkinson's above what simple baselines predict? (Trial entry is not efficacy.)

Positives   : Hetionet compounds with a phase 2/3 interventional trial for the disease starting on/after --cutoff
              (default 2017), that are not already Hetionet treats/palliates pairs.
Universe    : all compounds except those known pairs (rest are unlabeled negatives; compounds with no binding
              data are ranked last).
Ranking     : un-triaged consensus (mean of DWPC and permutation-z percentiles, binding metapaths) - the same
              definition as 09_consensus_candidates.py. Needs data/processed/zscores_w<w>.npz from 08.
Baselines   : degree (# bound genes); pre-cutoff trials of the compound for the same disease; optional
              pre-cutoff trials of the compound for ANY disease (--any-disease-counts, ~30 min of API calls);
              optional Rephetio probabilities (--rephetio TSV); random = AUROC 0.5.
Reported    : AUROC, hit@5%/10%, average precision, precision@20, paired bootstrap differences over positives,
              a novelty-only analysis (no pre-cutoff trial of that pair), sensitivity to cutoff year and phase,
              hub / non-therapeutic flags for the top 20, and the outcome class declared in the plan.

Outputs (outputs/tables/): ctgov_trials_<slug>.csv, ctgov_unmatched_<slug>.csv, temporal_positives_<slug>.csv,
    temporal_top20_<slug>.csv, temporal_summary.csv, temporal_sensitivity.csv
Raw API results are cached in data/raw/ctgov_<slug>.json (with the retrieval date) so reruns are offline.

Usage:
    python3 src/13_temporal_validation.py                       # downloads (cached) and evaluates
    python3 src/13_temporal_validation.py --offline             # cache only
    python3 src/13_temporal_validation.py --any-disease-counts  # adds the pre-cutoff any-disease baseline
Run `git tag frozen-v1` BEFORE the first run (see docs/validation_plan.md).
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import date

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.metrics import average_precision_score

from dwpc import MECHANISTIC, aggregate_vector, build_matrices, compute_dwpc, coverage, known_mask
from utils import DATA_PROC, DATA_RAW, OUT_TAB, ROOT, find_nodes, load_processed

API = "https://clinicaltrials.gov/api/v2/studies"
UA = "hetionet-repurposing-temporal-validation/1.0 (portfolio project)"
FIELDS = "NCTId,BriefTitle,OverallStatus,StartDate,Phase,StudyType,InterventionType,InterventionName,Condition"
DISEASES = {"alzheimer": ("Alzheimer Disease", "alzheimer"), "parkinson": ("Parkinson Disease", "parkinson")}
BIND = ["CbGaD", "CbGiGaD", "CbGpPWpGaD"]
INTERVENTION_TYPES = {"DRUG", "DIETARY_SUPPLEMENT", "BIOLOGICAL", "COMBINATION_PRODUCT"}
EXCLUDED_CATEGORIES = {"endogenous", "diagnostic", "cytotoxic_or_metal", "risk_factor", "wrong_direction_likely"}
STOP = {"hydrochloride", "hcl", "sodium", "potassium", "calcium", "mesylate", "tartrate", "maleate", "sulfate",
        "sulphate", "acetate", "citrate", "phosphate", "succinate", "besylate", "fumarate", "bromide", "chloride",
        "hydrobromide", "tablet", "tablets", "capsule", "capsules", "patch", "transdermal", "oral", "injection",
        "infusion", "extended", "release", "er", "xr", "sr", "mg", "dose", "low", "high", "placebo", "solution",
        "cream", "gel", "spray", "intranasal", "nasal", "iv", "intravenous", "subcutaneous", "daily"}


# --------------------------------------------------------------------------- ClinicalTrials.gov
def http_json(url: str, retries: int = 4) -> dict:
    for k in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503, 504) and k < retries - 1:
                time.sleep(2 ** (k + 1))
                continue
            raise
        except (urllib.error.URLError, TimeoutError):
            if k < retries - 1:
                time.sleep(2 ** (k + 1))
                continue
            raise
    raise RuntimeError("unreachable")


def fetch_all(api: str, cond: str, sleep: float = 0.4) -> list[dict]:
    """All studies for a condition query (client-side filtering later). Retries without `fields` on HTTP 400."""
    use_fields = True
    studies, token = [], None
    while True:
        params = {"query.cond": cond, "pageSize": "1000", "format": "json"}
        if use_fields:
            params["fields"] = FIELDS
        if token:
            params["pageToken"] = token
        try:
            data = http_json(api + "?" + urllib.parse.urlencode(params))
        except urllib.error.HTTPError as e:
            if e.code == 400 and use_fields:
                print("  HTTP 400 with `fields`; retrying without it (larger download)")
                use_fields, studies, token = False, [], None
                continue
            raise
        studies += data.get("studies", [])
        token = data.get("nextPageToken")
        print(f"  fetched {len(studies):,} studies", end="\r", flush=True)
        if not token:
            break
        time.sleep(sleep)
    print()
    return studies


def load_or_fetch(slug: str, cond: str, api: str, refresh: bool, offline: bool) -> tuple[list[dict], str]:
    path = DATA_RAW / f"ctgov_{slug}.json"
    if path.exists() and not refresh:
        blob = json.loads(path.read_text())
        return blob["studies"], blob["retrieved"]
    if offline:
        raise SystemExit(f"--offline but {path} is missing")
    print(f"Downloading ClinicalTrials.gov studies for '{cond}' ...")
    studies = fetch_all(api, cond)
    blob = {"retrieved": date.today().isoformat(), "condition": cond, "api": api, "n": len(studies),
            "studies": studies}
    path.write_text(json.dumps(blob))
    return studies, blob["retrieved"]


def parse_studies(studies: list[dict]) -> pd.DataFrame:
    rows = []
    for st in studies:
        ps = st.get("protocolSection", {})
        start = (ps.get("statusModule", {}).get("startDateStruct") or {}).get("date")
        year = int(start[:4]) if start and re.match(r"^\d{4}", start) else None
        ivs = [(i.get("type", ""), i.get("name", "")) for i in
               (ps.get("armsInterventionsModule", {}).get("interventions") or [])]
        rows.append({"nct": ps.get("identificationModule", {}).get("nctId"),
                     "title": ps.get("identificationModule", {}).get("briefTitle", ""),
                     "status": ps.get("statusModule", {}).get("overallStatus", ""),
                     "year": year, "study_type": ps.get("designModule", {}).get("studyType", ""),
                     "phases": ps.get("designModule", {}).get("phases") or [], "interventions": ivs})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- name matching
def norm(s: str) -> str:
    s = str(s).lower().replace("'", "")
    s = re.sub(r"\(.*?\)", " ", s)
    s = re.sub(r"[^a-z0-9\- ]+", " ", s)
    toks = [t for t in s.split() if t and t not in STOP and not re.fullmatch(r"\d+(\.\d+)?(mg|mcg|g|ml|%)?", t)]
    return " ".join(toks)


def build_index(comp_names: list[str], alias_path) -> tuple[dict, list[str]]:
    idx = {}
    for i, n in enumerate(comp_names):
        k = norm(n)
        if k and k not in idx:
            idx[k] = i
    missing = []
    if alias_path.exists():
        lower = {n.lower(): i for i, n in enumerate(comp_names)}
        for r in pd.read_csv(alias_path).itertuples():
            t = str(r.hetionet_name).lower()
            if t in lower:
                idx[norm(r.alias)] = lower[t]
            else:
                missing.append(f"{r.alias} -> {r.hetionet_name}")
    return idx, missing


def match_name(name: str, idx: dict) -> set[int]:
    if norm(name) in idx:
        return {idx[norm(name)]}
    parts = re.split(r"\s*(?:\+|/|,|;|\band\b|\bplus\b|\bwith\b)\s*", name.lower())
    return {idx[norm(p)] for p in parts if norm(p) in idx}


def trial_table(trials: pd.DataFrame, idx: dict) -> tuple[pd.DataFrame, pd.Series]:
    rows, unmatched = [], []
    for t in trials.itertuples():
        if t.study_type != "INTERVENTIONAL" or t.year is None:
            continue
        for typ, name in t.interventions:
            if typ not in INTERVENTION_TYPES or "placebo" in name.lower():
                continue
            m = match_name(name, idx)
            if not m:
                unmatched.append(name)
            for c in m:
                rows.append((t.nct, t.year, tuple(t.phases), c, name))
    tt = pd.DataFrame(rows, columns=["nct", "year", "phases", "c", "name"]).drop_duplicates(["nct", "c"])
    return tt, pd.Series(unmatched, dtype=str)


# --------------------------------------------------------------------------- scoring / metrics
def pct(v: np.ndarray) -> np.ndarray:
    out = np.full(len(v), np.nan)
    ok = np.isfinite(v)
    out[ok] = rankdata(v[ok], method="average") / ok.sum()
    return out


def consensus_vector(obs, Z, cov, d, min_nonzero):
    s_d = aggregate_vector(obs, cov, d, BIND, "mean", None, min_nonzero)[0]
    s_z = aggregate_vector(Z, cov, d, BIND, "mean", None, 0, signed_all=True)[0]
    return (pct(s_d) + pct(s_z)) / 2


def method_stats(v: np.ndarray, pos: np.ndarray, uni: np.ndarray):
    vu, pu = v[uni].astype(float), pos[uni]
    if pu.sum() == 0 or (~pu).sum() == 0 or np.isnan(vu).all():
        return None
    vu = np.where(np.isnan(vu), np.nanmin(vu) - 1.0, vu)
    neg = np.sort(vu[~pu])
    s = vu[pu]
    lo, hi = np.searchsorted(neg, s, "left"), np.searchsorted(neg, s, "right")
    u = (lo + 0.5 * (hi - lo)) / len(neg)
    rd = rankdata(-vu, method="average") / len(vu)
    top = np.argsort(-vu, kind="stable")[:20]
    return {"u": u, "auroc": float(u.mean()), "hit5": float((rd[pu] <= 0.05).mean()),
            "hit10": float((rd[pu] <= 0.10).mean()), "ap": float(average_precision_score(pu, vu)),
            "p20": float(pu[top].mean())}


def boot_diff(u_a: np.ndarray, u_b: np.ndarray, n_boot: int, rng) -> tuple[float, float, float]:
    d = u_a - u_b
    if len(d) < 2:
        return float(d.mean()), np.nan, np.nan
    m = d[rng.integers(0, len(d), (n_boot, len(d)))].mean(axis=1)
    return float(d.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def git_freeze_check() -> None:
    try:
        tags = subprocess.run(["git", "tag", "--list", "frozen-v1"], cwd=ROOT, capture_output=True, text=True,
                              timeout=10).stdout.split()
        if not tags:
            print("WARNING: git tag 'frozen-v1' not found. Create it before the first run "
                  "(docs/validation_plan.md) so changes after seeing results are visible.")
            return
        diff = subprocess.run(["git", "diff", "--name-only", "frozen-v1", "--", "src", "data/triage.csv"],
                              cwd=ROOT, capture_output=True, text=True, timeout=10).stdout.split()
        diff = [f for f in diff if not f.endswith("13_temporal_validation.py")]
        if diff:
            print("WARNING: files changed since frozen-v1 (report any effect as exploratory):", ", ".join(diff))
        else:
            print("Freeze check: pipeline unchanged since frozen-v1.")
    except Exception:
        print("(git not available; freeze check skipped)")


def load_rephetio(path: str, comp_ids: list[str], disease_id: str):
    df = pd.read_csv(path, sep="\t")
    low = {c.lower(): c for c in df.columns}
    ccol = next((low[c] for c in low if "drugbank" in c), None) or \
        next((low[c] for c in low if "compound" in c and "id" in c), None)
    dcol = next((low[c] for c in low if "doid" in c), None) or next((low[c] for c in low if "disease" in c and "id" in c), None)
    pcol = next((low[c] for c in low if "prob" in c), None)
    if not (ccol and dcol and pcol):
        print(f"  Rephetio file: could not identify columns in {list(df.columns)}; skipping comparator")
        return None
    doid = disease_id.split("::")[-1]
    sub = df[df[dcol].astype(str).str.contains(doid, regex=False)]
    prob = {str(r[ccol]).split("::")[-1]: r[pcol] for _, r in sub.iterrows()}
    return np.array([prob.get(i.split("::")[-1], np.nan) for i in comp_ids], dtype=float)


def fetch_any_counts(api: str, names: list[str], cutoff: int, sleep: float = 1.3):
    path = DATA_RAW / f"ctgov_anycounts_before{cutoff}.json"
    cache = json.loads(path.read_text()) if path.exists() else {}
    todo = [n for n in names if n not in cache]
    if todo:
        print(f"Any-disease counts: {len(todo)} compounds to query (~{len(todo) * sleep / 60:.0f} min)")
    for k, n in enumerate(todo, 1):
        params = {"query.intr": n, "filter.advanced": f"AREA[StartDate]RANGE[MIN,{cutoff - 1}-12-31]",
                  "countTotal": "true", "pageSize": "1", "format": "json", "fields": "NCTId"}
        try:
            cache[n] = int(http_json(api + "?" + urllib.parse.urlencode(params)).get("totalCount", 0))
        except urllib.error.HTTPError as e:
            print(f"  API error {e.code} for '{n}'; any-disease baseline aborted")
            return None
        time.sleep(sleep)
        if k % 25 == 0:
            path.write_text(json.dumps(cache))
            print(f"  {k}/{len(todo)}", end="\r", flush=True)
    path.write_text(json.dumps(cache))
    return cache


# --------------------------------------------------------------------------- main
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--w", type=float, default=0.6)
    ap.add_argument("--cutoff", type=int, default=2017)
    ap.add_argument("--min-nonzero", type=int, default=20)
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--api", default=API)
    ap.add_argument("--refresh", action="store_true")
    ap.add_argument("--offline", action="store_true")
    ap.add_argument("--any-disease-counts", action="store_true")
    ap.add_argument("--count-sleep", type=float, default=1.3, help="seconds between count queries (API allows ~50/min)")
    ap.add_argument("--rephetio", default=None, help="TSV of Rephetio probabilities (best-effort column detection)")
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    git_freeze_check()
    zpath = DATA_PROC / f"zscores_w{args.w:g}.npz"
    if not zpath.exists():
        raise SystemExit(f"{zpath} not found - run 08_permutation_null.py first")
    Z = {k: np.load(zpath)[k] for k in MECHANISTIC}

    nodes, edges = load_processed()
    mats = build_matrices(nodes, edges)
    A = mats.A
    cov = coverage(A)
    comp = mats.tables["Compound"]
    comp_names, comp_ids = comp["name"].tolist(), comp["id"].tolist()
    deg = np.asarray(A["CbG"].sum(axis=1)).ravel().astype(float)
    obs = compute_dwpc(A, args.w)
    idx, alias_missing = build_index(comp_names, ROOT / "data" / "drug_synonyms.csv")
    if alias_missing:
        print("Synonym targets not found in Hetionet (ignored):", "; ".join(alias_missing))

    triage = {}
    tp = ROOT / "data" / "triage.csv"
    if tp.exists():
        triage = dict(zip(pd.read_csv(tp)["compound"].str.lower(), pd.read_csv(tp)["category"]))

    any_counts = None
    if args.any_disease_counts:
        any_counts = fetch_any_counts(args.api, [n for n, c in zip(comp_names, cov["CbGaD"]) if c], args.cutoff, args.count_sleep)

    data = {}
    for slug, (cond, hname) in DISEASES.items():
        hit = find_nodes(nodes, hname, kind="Disease").iloc[0]
        d = mats.local_index(int(hit["idx"]))
        studies, retrieved = load_or_fetch(slug, cond, args.api, args.refresh, args.offline)
        trials = parse_studies(studies)
        tt, unmatched = trial_table(trials, idx)
        n_int = int((trials["study_type"] == "INTERVENTIONAL").sum())
        print(f"{hit['name']}: {len(trials):,} studies (retrieved {retrieved}), {n_int:,} interventional, "
              f"{tt['c'].nunique()} Hetionet compounds matched; unmatched intervention mentions: {len(unmatched):,}")
        unmatched.value_counts().head(60).rename_axis("intervention").reset_index(name="n").to_csv(
            OUT_TAB / f"ctgov_unmatched_{slug}.csv", index=False)
        tt.assign(compound=[comp_names[c] for c in tt["c"]]).to_csv(OUT_TAB / f"ctgov_trials_{slug}.csv", index=False)
        vec = {"consensus": consensus_vector(obs, Z, cov, d, args.min_nonzero), "degree": deg}
        if args.rephetio:
            rv = load_rephetio(args.rephetio, comp_ids, hit["id"])
            if rv is not None:
                vec["rephetio"] = rv
        data[slug] = dict(d=d, name=hit["name"], tt=tt, vec=vec, known=known_mask(A, d))

    def config(cutoff: int, phases: set, novelty: bool = False):
        pooled, per = {}, {}
        for slug, D in data.items():
            tt = D["tt"]
            pre = tt[tt["year"] < cutoff].groupby("c")["nct"].nunique()
            post = tt[(tt["year"] >= cutoff) & tt["phases"].apply(lambda ph: any(p in phases for p in ph))]
            pos = np.zeros(len(comp_names), bool)
            pos[post["c"].unique()] = True
            pos &= ~D["known"]
            pre_vec = np.zeros(len(comp_names))
            pre_vec[pre.index.to_numpy()] = pre.to_numpy()
            uni = ~D["known"]
            if novelty:
                uni = uni & (pre_vec == 0)
                pos &= uni
            vecs = dict(D["vec"])
            vecs["trials_disease"] = pre_vec
            if any_counts is not None and cutoff == args.cutoff:
                vecs["trials_any"] = np.array([any_counts.get(n, np.nan) for n in comp_names], dtype=float)
            stats = {m: method_stats(v, pos, uni) for m, v in vecs.items()}
            per[slug] = dict(stats=stats, pos=pos, uni=uni, pre=pre_vec, post=post)
            for m, s in stats.items():
                if s is not None:
                    pooled.setdefault(m, []).append(s["u"])
        return {m: np.concatenate(v) for m, v in pooled.items()}, per

    def report(cutoff, phases, novelty=False, label=""):
        pooled, per = config(cutoff, phases, novelty)
        n_pos = len(pooled.get("consensus", []))
        rows = []
        for scope, getter in [("pooled", None)] + [(s, s) for s in data]:
            for m in ["consensus", "degree", "trials_disease", "trials_any", "rephetio"]:
                if getter is None:
                    u = pooled.get(m)
                    st = [per[s]["stats"].get(m) for s in data if per[s]["stats"].get(m) is not None]
                    if u is None:
                        continue
                    row = dict(scope=scope, method=m, n_pos=len(u), auroc=float(u.mean()),
                               hit5=np.mean(np.concatenate([[x["hit5"]] * len(x["u"]) for x in st])),
                               hit10=np.mean(np.concatenate([[x["hit10"]] * len(x["u"]) for x in st])),
                               ap=np.mean([x["ap"] for x in st]), p20=np.mean([x["p20"] for x in st]))
                else:
                    s = per[getter]["stats"].get(m)
                    if s is None:
                        continue
                    row = dict(scope=scope, method=m, n_pos=len(s["u"]), auroc=s["auroc"], hit5=s["hit5"],
                               hit10=s["hit10"], ap=s["ap"], p20=s["p20"])
                rows.append(row)
        diffs = {}
        for b in ("degree", "trials_disease", "trials_any", "rephetio"):
            if b in pooled and "consensus" in pooled and len(pooled[b]) == n_pos:
                diffs[b] = boot_diff(pooled["consensus"], pooled[b], args.boot, rng)
        return rows, diffs, per, n_pos

    phases_main = {"PHASE2", "PHASE3"}
    rows, diffs, per, n_pos = report(args.cutoff, phases_main)
    summ = pd.DataFrame(rows)
    print(f"\n=== Primary analysis: trials starting >= {args.cutoff}, phase 2/3 (pooled positives = {n_pos}) ===")
    for slug, D in data.items():
        print(f"  {D['name']}: {int(per[slug]['pos'].sum())} positives among {int(per[slug]['uni'].sum())} compounds "
              f"({int(D['known'].sum())} known treats/palliates removed)")
    print(summ.round(3).to_string(index=False))
    summ.to_csv(OUT_TAB / "temporal_summary.csv", index=False)
    print("\nPaired difference in AUROC, consensus minus baseline (95% bootstrap CI over positives):")
    for b, (dm, lo, hi) in diffs.items():
        print(f"  vs {b:<15} {dm:+.3f}  [{lo:+.3f}, {hi:+.3f}]")
    if "trials_any" not in diffs:
        print("  (any-disease trial-count baseline not run; use --any-disease-counts. Report this deviation.)")

    # novelty-only
    rows_n, diffs_n, per_n, n_pos_n = report(args.cutoff, phases_main, novelty=True)
    cn = next((r for r in rows_n if r["scope"] == "pooled" and r["method"] == "consensus"), None)
    print(f"\n=== Novelty-only (no pre-{args.cutoff} trial of the pair; pooled positives = {n_pos_n}) ===")
    if cn:
        print(f"  consensus AUROC {cn['auroc']:.3f}; " + "; ".join(
            f"vs {b} {dm:+.3f} [{lo:+.3f}, {hi:+.3f}]" for b, (dm, lo, hi) in diffs_n.items() if b != "trials_disease"))

    # sensitivity
    sens = []
    for cut in (2017, 2018, 2019):
        for ph, lab in ((phases_main, "phase 2/3"), ({"PHASE3"}, "phase 3")):
            _, df_, _, n_ = report(cut, ph)
            for b, (dm, lo, hi) in df_.items():
                sens.append(dict(cutoff=cut, phases=lab, n_pos=n_, baseline=b, delta=dm, ci_lo=lo, ci_hi=hi))
    sens = pd.DataFrame(sens)
    sens.to_csv(OUT_TAB / "temporal_sensitivity.csv", index=False)
    print("\n=== Sensitivity (pooled AUROC difference consensus - baseline) ===")
    print(sens.round(3).to_string(index=False))

    # top-20 diagnostics and tables
    hub_flag = cat_flag = False
    for slug, D in data.items():
        P = per[slug]
        cons = D["vec"]["consensus"]
        order = [i for i in np.argsort(-np.nan_to_num(cons, nan=-1), kind="stable")
                 if P["uni"][i] and P["pre"][i] == 0][:20]
        top = pd.DataFrame({"compound": [comp_names[i] for i in order], "consensus": cons[order],
                            "n_bound_genes": deg[order], "became_trial_positive": P["pos"][order],
                            "category": [triage.get(comp_names[i].lower(), "") for i in order]})
        top.to_csv(OUT_TAB / f"temporal_top20_{slug}.csv", index=False)
        med_top, med_all = np.median(deg[order]), np.median(deg[cov["CbGaD"]])
        share = float(top["category"].isin(EXCLUDED_CATEGORIES).mean())
        hub_flag |= med_top > 3 * med_all
        cat_flag |= share > 0.25
        print(f"\n{D['name']}: top-20 novel candidates: {int(top['became_trial_positive'].sum())} became phase 2/3 "
              f"trial positives; median genes bound {med_top:.0f} vs {med_all:.0f}; "
              f"{share:.0%} in non-therapeutic triage categories")
        post = P["post"]
        pos_idx = np.where(P["pos"])[0]
        rank_in_uni = {}
        vu = np.where(np.isnan(cons), np.nanmin(cons) - 1, cons)
        rk = rankdata(-vu[P["uni"]], method="min")
        rank_in_uni = dict(zip(np.where(P["uni"])[0], rk))
        pt = pd.DataFrame({
            "compound": [comp_names[i] for i in pos_idx],
            "consensus_rank": [int(rank_in_uni[i]) for i in pos_idx],
            "universe_size": int(P["uni"].sum()),
            "n_bound_genes": deg[pos_idx], "pre_trials_same_disease": P["pre"][pos_idx],
            "n_post_trials": [int(post[post["c"] == i]["nct"].nunique()) for i in pos_idx],
            "example_ncts": [", ".join(post[post["c"] == i]["nct"].head(3)) for i in pos_idx],
            "matched_from": ["; ".join(sorted(set(post[post["c"] == i]["name"]))[:3]) for i in pos_idx]})
        pt.sort_values("consensus_rank").to_csv(OUT_TAB / f"temporal_positives_{slug}.csv", index=False)

    # outcome class (as declared in docs/validation_plan.md)
    print("\n=== Outcome class ===")
    if n_pos < 20:
        outcome = "INCONCLUSIVE - fewer than 20 pooled positives"
    else:
        beats = {b: lo > 0 for b, (dm, lo, hi) in diffs.items()}
        nov_ok = bool(diffs_n.get("degree") and diffs_n["degree"][1] > 0) if n_pos_n >= 10 else None
        if not beats.get("degree", False):
            outcome = "NULL - does not beat the degree baseline"
        elif not all(beats.values()):
            outcome = "PARTIAL - beats degree but not every baseline (graph tracks existing research interest)"
        elif nov_ok is False:
            outcome = "LEAKAGE SUSPECT - beats all baselines overall but not in the novelty-only analysis"
        else:
            outcome = "SUCCESS - beats every available baseline" + ("" if nov_ok else " (novelty-only not assessable)")
    print(outcome)
    if hub_flag:
        print("FLAG: hub failure - top-20 candidates bind more than 3x the typical number of genes")
    if cat_flag:
        print("FLAG: more than 25% of top-20 candidates fall in non-therapeutic triage categories")
    print("Direction of effect (contraindications) is not assessed here; see the plan's optional add-on.")


if __name__ == "__main__":
    main()
