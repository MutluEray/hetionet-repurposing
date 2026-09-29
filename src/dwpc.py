"""Sparse-matrix DWPC engine for Hetionet (compound x disease scoring).

DWPC = degree-weighted path count (Himmelstein et al. 2017, Project Rephetio).
Each relation's adjacency matrix is degree-weighted as  D_row^-w * A * D_col^-w
(degrees taken with respect to that relation), and a metapath score matrix is the
product of the weighted matrices along the metapath. w = 0 gives raw path counts.

Metapaths (compound -> disease):
    CbGaD        compound binds gene associated with disease
    CbGiGaD      compound binds gene that interacts with a disease-associated gene
    CbGpPWpGaD   compound binds gene sharing a pathway with a disease-associated gene
                 (degenerate paths that revisit the same gene are subtracted)
    SigNet       expression-signature reversal minus concordance:
                 (CuGdD + CdGuD) - (CuGuD + CdGdD)
    CrCtD        compound resembles a compound that treats the disease
    CtDrD        compound treats a disease that resembles the target disease
Evidence groups: binding (CbGaD, CbGiGaD, CbGpPWpGaD), signature (SigNet),
similarity (CrCtD, CtDrD).

Aggregation (v2):
  * Per metapath, scores become percentile ranks among the compounds that have input
    data for that metapath (coverage-aware). Non-negative metapaths use "min" ranking so
    zero-evidence compounds get percentile 0 (not a neutral ~0.5); the signed SigNet keeps
    average ranking so zero sits mid-scale.
  * A metapath is dropped for a disease if fewer than `min_nonzero` covered compounds have
    non-zero evidence (near-empty metapaths only inject arbitrary offsets).
  * mode="mean":   mean percentile over the metapaths that remain.
    mode="groups": mean within each evidence group first, then across groups (prevents the
                   three highly correlated binding-based metapaths from triple-counting).
"""
from __future__ import annotations

import re
from dataclasses import dataclass

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import rankdata

METAPATHS = ["CbGaD", "CbGiGaD", "CbGpPWpGaD", "SigNet", "CrCtD", "CtDrD"]
MECHANISTIC = ["CbGaD", "CbGiGaD", "CbGpPWpGaD", "SigNet"]
SIGNED = {"SigNet"}
GROUPS = {
    "binding": ["CbGaD", "CbGiGaD", "CbGpPWpGaD"],
    "signature": ["SigNet"],
    "similarity": ["CrCtD", "CtDrD"],
}
EPS = 1e-12

REL = {  # metaedge: (source kind, target kind)
    "CbG": ("Compound", "Gene"), "CuG": ("Compound", "Gene"), "CdG": ("Compound", "Gene"),
    "DaG": ("Disease", "Gene"), "DuG": ("Disease", "Gene"), "DdG": ("Disease", "Gene"),
    "GpPW": ("Gene", "Pathway"), "GiG": ("Gene", "Gene"),
    "CrC": ("Compound", "Compound"), "CtD": ("Compound", "Disease"),
    "CpD": ("Compound", "Disease"), "DrD": ("Disease", "Disease"),
}
SYMMETRIC = {"GiG", "CrC", "DrD"}  # stored once per pair in Hetionet -> symmetrize


@dataclass
class HetMats:
    A: dict                 # metaedge -> scipy CSR (binary; symmetric relations symmetrized)
    tables: dict            # kind -> DataFrame(id, name, idx) in matrix order
    loc: np.ndarray         # global node idx -> row/col position within its own kind
    dropped_pathways: list | None = None   # names removed by drop_disease_pathways

    def local_index(self, global_idx: int) -> int:
        return int(self.loc[global_idx])


def _norm(s: str) -> str:
    s = str(s).lower().replace("'", "").replace("\u2019", "")
    return re.sub(r"[^a-z0-9]+", " ", s).strip()


def disease_pathway_mask(tables: dict, min_len: int = 5) -> np.ndarray:
    """Boolean mask over pathways whose normalised name contains a Hetionet disease name
    (e.g. 'Alzheimers Disease', 'Parkinsons Disease Pathway', 'Bladder Cancer').
    Such pathways make the pathway metapath partly circular for the disease they are named after."""
    dnames = [n for n in (_norm(x) for x in tables["Disease"]["name"]) if len(n) >= min_len]
    return np.array([any(dn in _norm(p) for dn in dnames) for p in tables["Pathway"]["name"]])


def build_matrices(nodes: pd.DataFrame, edges: pd.DataFrame,
                   drop_disease_pathways: bool = False) -> HetMats:
    kind_arr = nodes["kind"].astype(str).to_numpy()
    loc = np.full(len(nodes), -1, dtype=np.int64)
    tables = {}
    for kind in sorted({k for pair in REL.values() for k in pair}):
        sub = nodes.loc[kind_arr == kind, ["id", "name", "idx"]].reset_index(drop=True)
        loc[sub["idx"].to_numpy()] = np.arange(len(sub))
        tables[kind] = sub
    me = edges["metaedge"].astype(str).to_numpy()
    src = edges["src_idx"].to_numpy()
    dst = edges["dst_idx"].to_numpy()
    A = {}
    for m, (sk, dk) in REL.items():
        k = me == m
        M = sp.csr_matrix((np.ones(int(k.sum())), (loc[src[k]], loc[dst[k]])),
                          shape=(len(tables[sk]), len(tables[dk])))
        if m in SYMMETRIC:
            M = M + M.T
        M = M.tocsr()
        M.data[:] = 1.0
        A[m] = M
    dropped = None
    if drop_disease_pathways:
        mask = disease_pathway_mask(tables)
        M = (A["GpPW"] @ sp.diags((~mask).astype(float))).tocsr()
        M.eliminate_zeros()
        A["GpPW"] = M
        dropped = tables["Pathway"]["name"][mask].tolist()
    return HetMats(A=A, tables=tables, loc=loc, dropped_pathways=dropped)


def _dw(M: sp.spmatrix, w: float) -> sp.csr_matrix:
    if w == 0:
        return M.astype(np.float64).tocsr()
    dr = np.asarray(M.sum(axis=1)).ravel()
    dc = np.asarray(M.sum(axis=0)).ravel()
    ir = np.zeros_like(dr, dtype=np.float64)
    ic = np.zeros_like(dc, dtype=np.float64)
    ir[dr > 0] = dr[dr > 0] ** (-w)
    ic[dc > 0] = dc[dc > 0] ** (-w)
    return (sp.diags(ir) @ M @ sp.diags(ic)).tocsr()


def _chain(*mats) -> np.ndarray:
    """Right-to-left product; the last factor is made dense so every step stays small."""
    last = mats[-1]
    acc = last.toarray() if sp.issparse(last) else np.asarray(last)
    for m in reversed(mats[:-1]):
        acc = m @ acc
    return np.asarray(acc)


def compute_dwpc(A: dict, w: float) -> dict:
    """Return {metapath: dense array (n_compounds x n_diseases)} for damping exponent w."""
    CbG, CuG, CdG = (_dw(A[k], w) for k in ("CbG", "CuG", "CdG"))
    DaG, DuG, DdG = (_dw(A[k], w) for k in ("DaG", "DuG", "DdG"))
    GiG, GP = _dw(A["GiG"], w), _dw(A["GpPW"], w)
    CrC, CtD, DrD = _dw(A["CrC"], w), _dw(A["CtD"], w), _dw(A["DrD"], w)

    out = {}
    out["CbGaD"] = _chain(CbG, DaG.T)
    out["CbGiGaD"] = _chain(CbG, GiG, DaG.T)

    pw = _chain(CbG, GP, GP.T, DaG.T)
    s = np.asarray(GP.multiply(GP).sum(axis=1)).ravel()      # weight of g -> p -> g (same gene)
    pw = pw - _chain(CbG, sp.diags(s), DaG.T)                 # remove non-simple paths
    out["CbGpPWpGaD"] = np.clip(pw, 0, None)

    rev = _chain(CuG, DdG.T) + _chain(CdG, DuG.T)
    conc = _chain(CuG, DuG.T) + _chain(CdG, DdG.T)
    out["SigNet"] = rev - conc
    out["SigRev"], out["SigConc"] = rev, conc                 # kept for explanation

    out["CrCtD"] = _chain(CrC, CtD)
    out["CtDrD"] = _chain(CtD, DrD)

    for k, v in out.items():                                  # remove floating-point dust
        v = np.array(v, dtype=np.float64)
        v[np.abs(v) < EPS] = 0.0
        out[k] = v
    return out


def coverage(A: dict) -> dict:
    """Which compounds have any data feeding each metapath (bool arrays over compounds)."""
    has = lambda k: np.asarray(A[k].sum(axis=1)).ravel() > 0
    return {
        "CbGaD": has("CbG"), "CbGiGaD": has("CbG"), "CbGpPWpGaD": has("CbG"),
        "SigNet": has("CuG") | has("CdG"), "CrCtD": has("CrC"), "CtDrD": has("CtD"),
    }


def known_mask(A: dict, d: int) -> np.ndarray:
    """Compounds with a treats or palliates edge to disease column d."""
    return (A["CtD"][:, d].toarray().ravel() > 0) | (A["CpD"][:, d].toarray().ravel() > 0)


def aggregate_vector(scores: dict, cov: dict, d: int, mps: list, mode: str = "mean",
                     min_covered: int | None = None, min_nonzero: int = 20,
                     signed_all: bool = False):
    """Aggregate metapath scores for disease column d.

    signed_all=True treats every metapath as signed (average-rank percentiles, non-zero
    filter on != 0); use it for z-scores against a permutation null.

    Returns (score vector over compounds [NaN = unranked], percentile matrix
    (compounds x len(mps), NaN where uncovered/dropped), list of metapaths actually used).
    """
    C = scores[mps[0]].shape[0]
    pcts = np.full((C, len(mps)), np.nan)
    used = []
    for k, mp in enumerate(mps):
        m = cov[mp]
        if not m.any():
            continue
        s = scores[mp][m, d].astype(float)
        signed = signed_all or mp in SIGNED
        n_nz = int((s != 0).sum()) if signed else int((s > 0).sum())
        if n_nz < min_nonzero:
            continue
        if signed:
            p = rankdata(s, method="average") / len(s)
        else:
            p = (rankdata(s, method="min") - 1.0) / len(s)
        pcts[m, k] = p
        used.append(mp)

    if min_covered is None:
        min_covered = 2 if mode == "mean" else 1
    if mode == "mean":
        cnt = np.isfinite(pcts).sum(axis=1)
        score = np.nansum(pcts, axis=1) / np.maximum(cnt, 1)
    elif mode == "groups":
        cols_out = []
        for members in GROUPS.values():
            cols = [mps.index(x) for x in members if x in mps]
            if not cols:
                continue
            block = pcts[:, cols]
            c = np.isfinite(block).sum(axis=1)
            cols_out.append(np.where(c > 0, np.nansum(block, axis=1) / np.maximum(c, 1), np.nan))
        G = np.column_stack(cols_out) if cols_out else np.full((C, 1), np.nan)
        cnt = np.isfinite(G).sum(axis=1)
        score = np.nansum(G, axis=1) / np.maximum(cnt, 1)
    else:
        raise ValueError("mode must be 'mean' or 'groups'")
    score = np.where(cnt >= min_covered, score, np.nan)
    return score, pcts, used


def rank_disease(mats: HetMats, scores: dict, cov: dict, d: int, mps: list,
                 min_covered: int | None = None, mode: str = "mean", min_nonzero: int = 20,
                 signed_all: bool = False):
    """Rank all compounds for disease column d. Returns (DataFrame, metapaths used)."""
    score, pcts, used = aggregate_vector(scores, cov, d, mps, mode, min_covered, min_nonzero,
                                         signed_all=signed_all)
    comp = mats.tables["Compound"]
    df = pd.DataFrame({"compound_id": comp["id"], "compound": comp["name"],
                       "cidx": np.arange(len(comp)), "score": score,
                       "n_covered": np.isfinite(pcts).sum(axis=1),
                       "known": known_mask(mats.A, d)})
    for k, mp in enumerate(mps):
        df[f"pct_{mp}"] = pcts[:, k]
    for mp in mps:
        df[f"dwpc_{mp}"] = scores[mp][:, d]
    df["rank"] = df["score"].rank(ascending=False, method="min")
    df = df.sort_values("score", ascending=False, na_position="last").reset_index(drop=True)
    return df, used
