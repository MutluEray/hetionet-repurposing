"""Evaluation helpers shared by 07_evaluate.py and 08_permutation_null.py."""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr, wilcoxon

warnings.filterwarnings("ignore")


def metrics(v: np.ndarray, pos: np.ndarray):
    """AUROC, median rank percentile of positives, hit@5%, hit@10%. NaN scores rank last."""
    if np.isnan(v).all():
        return np.nan, np.nan, np.nan, np.nan
    v = np.where(np.isnan(v), np.nanmin(v) - 1.0, v)
    n, npos = len(v), int(pos.sum())
    nneg = n - npos
    r = rankdata(v)
    auc = (r[pos].sum() - npos * (npos + 1) / 2) / (npos * nneg)
    rd = rankdata(-v, method="average") / n
    return auc, float(np.median(rd[pos])), float((rd[pos] <= 0.05).mean()), float((rd[pos] <= 0.10).mean())


def rho(v: np.ndarray, deg: np.ndarray) -> float:
    """Spearman correlation between a score vector and compound degree (finite entries only)."""
    ok = np.isfinite(v)
    if ok.sum() < 5 or np.ptp(v[ok]) == 0 or np.ptp(deg[ok]) == 0:
        return np.nan
    return float(spearmanr(v[ok], deg[ok])[0])


def summarize(per: pd.DataFrame, n_boot: int, rng: np.random.Generator,
              base: str = "degree_CbG") -> pd.DataFrame:
    """Macro-average per-disease results: mean AUROC + bootstrap CI over diseases,
    win-rate and paired Wilcoxon p-value versus the baseline method."""
    piv = per.pivot(index="disease_id", columns="method", values="auroc")
    n_unranked = piv.isna().sum()
    X = piv.fillna(0.5).to_numpy()
    n = X.shape[0]
    means = X[rng.integers(0, n, (n_boot, n))].mean(axis=1)
    lo, hi = np.percentile(means, [2.5, 97.5], axis=0)
    base_col = X[:, list(piv.columns).index(base)]
    agg = per.groupby("method").agg(mean_median_rank_pct=("median_rank_pct", "mean"),
                                    mean_hit5=("hit5", "mean"), rho_degree=("rho_degree", "mean"))
    out = []
    for j, m in enumerate(piv.columns):
        diff = X[:, j] - base_col
        try:
            p = wilcoxon(diff)[1] if np.any(diff != 0) else np.nan
        except Exception:
            p = np.nan
        out.append(dict(method=m, mean_auroc=X[:, j].mean(), ci_lo=lo[j], ci_hi=hi[j],
                        median_auroc=np.median(X[:, j]), wins_vs_degree=(diff > 0).mean(),
                        p_vs_degree=p, n_diseases=n, n_unranked=int(n_unranked[m]),
                        **agg.loc[m].to_dict()))
    return pd.DataFrame(out).sort_values("mean_auroc", ascending=False)
