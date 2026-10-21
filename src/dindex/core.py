"""Core D-index calculations.

The D-index contrasts two groups of taxa, a "plus" group P and a "minus" group M:

    D = (P - M) / (P + M) = 2 p - 1 = tanh(0.5 * ln(P / M)),   with p = P / (P + M)

D runs from -1 (only M detected) through 0 (equal abundance) to +1 (only P detected). It depends only on the ratio P/M, so it
is unchanged by sequencing depth, by normalization to relative abundance, or by any factor applied equally to both groups.
It is undefined when P + M = 0.
"""
from __future__ import annotations

from typing import Iterable, Sequence

import numpy as np
import pandas as pd
from scipy import stats

__all__ = ["d_index", "log_ratio", "d_from_log_ratio", "aggregate_groups", "d_from_table", "d_ci_counts", "d_ci_bootstrap"]


def d_index(plus, minus):
    """D = (plus - minus) / (plus + minus), element-wise. Returns NaN where plus + minus == 0.

    `plus` and `minus` may be scalars, arrays or pandas objects of non-negative abundances in the same units
    (reads, coverage, TPM or relative abundance).
    """
    p = np.asarray(plus, dtype=float)
    m = np.asarray(minus, dtype=float)
    if np.any(p < 0) or np.any(m < 0):
        raise ValueError("abundances must be non-negative")
    den = p + m
    with np.errstate(invalid="ignore", divide="ignore"):
        d = np.where(den > 0, (p - m) / np.where(den > 0, den, 1.0), np.nan)
    if isinstance(plus, pd.Series):
        return pd.Series(d, index=plus.index, name="D")
    return d if d.ndim else float(d)


def log_ratio(plus, minus, pseudocount: float | None = None):
    """ln(plus / minus). Zeros give +/-inf unless a pseudocount is added to both groups."""
    p = np.asarray(plus, dtype=float)
    m = np.asarray(minus, dtype=float)
    if np.any(p < 0) or np.any(m < 0):
        raise ValueError("abundances must be non-negative")
    if pseudocount is not None:
        if not np.isfinite(pseudocount) or pseudocount <= 0:
            raise ValueError("pseudocount must be a finite positive number")
        p, m = p + pseudocount, m + pseudocount
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.log(p) - np.log(m)


def d_from_log_ratio(lr):
    """Inverse link: D = tanh(lr / 2). Use to convert predictions of a log-ratio model back to the D scale."""
    return np.tanh(np.asarray(lr, dtype=float) / 2.0)


def aggregate_groups(abundance: pd.DataFrame, taxonomy: pd.Series | dict, plus: str | Sequence[str], minus: str | Sequence[str],
                     features_as: str = "columns") -> pd.DataFrame:
    """Sum feature abundances into the plus group, the minus group and the total over all mapped features.

    Parameters
    ----------
    abundance : samples x features table (``features_as="columns"``) or features x samples (``features_as="rows"``).
    taxonomy : mapping from feature ID to group label (e.g. order, genus, ecotype). Unmapped features count toward the total only
        if they appear in `taxonomy`; features absent from `taxonomy` are ignored and reported.
    plus, minus : a group label or a list of labels. The two sets must not overlap.

    Returns a samples table with columns ``plus``, ``minus``, ``total`` (all mapped features) and ``n_plus_features``,
    ``n_minus_features`` (features with non-zero abundance).
    """
    if features_as not in {"columns", "rows"}:
        raise ValueError("features_as must be 'columns' or 'rows'")
    tab = abundance if features_as == "columns" else abundance.T
    tax = pd.Series(taxonomy)
    plus_set, minus_set = _as_set(plus), _as_set(minus)
    if plus_set & minus_set:
        raise ValueError(f"plus and minus groups overlap: {sorted(plus_set & minus_set)}")
    mapped = [c for c in tab.columns if c in tax.index]
    if not mapped:
        raise ValueError("no abundance features found in taxonomy")
    groups = tax.reindex(mapped).astype(str)
    pcols = groups.index[groups.isin(plus_set)]
    mcols = groups.index[groups.isin(minus_set)]
    missing = sorted((plus_set | minus_set) - set(groups.unique()))
    if missing:
        raise ValueError(f"group labels not present in taxonomy for these features: {missing}")
    vals = tab[mapped].apply(pd.to_numeric, errors="raise")
    if (vals < 0).any().any():
        raise ValueError("abundances must be non-negative")
    out = pd.DataFrame({
        "plus": vals[pcols].sum(axis=1),
        "minus": vals[mcols].sum(axis=1),
        "total": vals.sum(axis=1),
        "n_plus_features": (vals[pcols] > 0).sum(axis=1),
        "n_minus_features": (vals[mcols] > 0).sum(axis=1),
    })
    out.attrs["unmapped_features"] = [c for c in tab.columns if c not in tax.index]
    return out


def d_from_table(abundance: pd.DataFrame, taxonomy: pd.Series | dict, plus, minus, features_as: str = "columns",
                 min_share: float | None = None, min_total: float | None = None, pseudocount: float | None = None) -> pd.DataFrame:
    """Compute D per sample from a feature table and a feature -> group mapping.

    Adds QC columns:
      ``share``  : (plus + minus) / total over all mapped features; low values mean D summarizes a minor part of the community.
      ``qc``     : "ok", "undefined" (plus + minus == 0), "low_share" (share < min_share) or "low_total" (plus + minus < min_total).
    ``log_ratio`` is ln(plus/minus), with an optional pseudocount; ``D`` never needs one.
    """
    if min_share is not None and not 0 <= min_share <= 1:
        raise ValueError("min_share must be between 0 and 1")
    if min_total is not None and min_total < 0:
        raise ValueError("min_total must be non-negative")
    agg = aggregate_groups(abundance, taxonomy, plus, minus, features_as)
    agg["D"] = d_index(agg["plus"], agg["minus"])
    agg["log_ratio"] = log_ratio(agg["plus"], agg["minus"], pseudocount)
    with np.errstate(invalid="ignore", divide="ignore"):
        agg["share"] = (agg["plus"] + agg["minus"]) / agg["total"].where(agg["total"] > 0)
    qc = np.full(len(agg), "ok", dtype=object)
    if min_total is not None:
        qc[(agg["plus"] + agg["minus"]).to_numpy() < min_total] = "low_total"
    if min_share is not None:
        qc[(agg["share"].fillna(0).to_numpy() < min_share) & (qc == "ok")] = "low_share"
    qc[((agg["plus"] + agg["minus"]) == 0).to_numpy()] = "undefined"
    agg["qc"] = qc
    return agg


def d_ci_counts(n_plus, n_minus, level: float = 0.95, prior: float = 0.5):
    """Credible interval for D from read counts assigned to each group.

    Uses p = P/(P+M) ~ Beta(n_plus + prior, n_minus + prior) (Jeffreys prior by default) and D = 2p - 1. Suitable when the
    inputs are counts of independent reads; for coverage or TPM values use :func:`d_ci_bootstrap`.
    Returns arrays (estimate, lower, upper); the estimate is the observed D.
    """
    p = np.asarray(n_plus, dtype=float)
    m = np.asarray(n_minus, dtype=float)
    if np.any(p < 0) or np.any(m < 0):
        raise ValueError("counts must be non-negative")
    if not 0 < level < 1:
        raise ValueError("level must be between 0 and 1")
    if not np.isfinite(prior) or prior <= 0:
        raise ValueError("prior must be a finite positive number")
    a = p + prior
    b = m + prior
    alpha = (1 - level) / 2
    lo = 2 * stats.beta.ppf(alpha, a, b) - 1
    hi = 2 * stats.beta.ppf(1 - alpha, a, b) - 1
    return d_index(n_plus, n_minus), lo, hi


def d_ci_bootstrap(abundance: pd.DataFrame, taxonomy, plus, minus, features_as: str = "columns", n_boot: int = 2000,
                   level: float = 0.95, seed: int | None = 0) -> pd.DataFrame:
    """Percentile interval for D by resampling features (e.g. phylotypes or genomes) with replacement within each group.

    This is a descriptive feature-composition sensitivity interval, not a sampling-error confidence interval. It measures
    how strongly D depends on the finite set of observed features within each group.
    """
    if features_as not in {"columns", "rows"}:
        raise ValueError("features_as must be 'columns' or 'rows'")
    if not isinstance(n_boot, (int, np.integer)) or n_boot < 2:
        raise ValueError("n_boot must be an integer of at least 2")
    if not 0 < level < 1:
        raise ValueError("level must be between 0 and 1")
    tab = abundance if features_as == "columns" else abundance.T
    tax = pd.Series(taxonomy).astype(str)
    pcols = [c for c in tab.columns if tax.get(c) in _as_set(plus)]
    mcols = [c for c in tab.columns if tax.get(c) in _as_set(minus)]
    if not pcols or not mcols:
        raise ValueError("taxonomy must contain at least one abundance feature from each group")
    P, M = tab[pcols].to_numpy(float), tab[mcols].to_numpy(float)
    if np.any(P < 0) or np.any(M < 0):
        raise ValueError("abundances must be non-negative")
    rng = np.random.default_rng(seed)
    boots = np.empty((n_boot, len(tab)))
    for i in range(n_boot):
        bp = P[:, rng.integers(0, P.shape[1], P.shape[1])].sum(1) if P.shape[1] else np.zeros(len(tab))
        bm = M[:, rng.integers(0, M.shape[1], M.shape[1])].sum(1) if M.shape[1] else np.zeros(len(tab))
        boots[i] = d_index(bp, bm)
    alpha = (1 - level) / 2
    with np.errstate(all="ignore"):
        return pd.DataFrame({"D": d_index(P.sum(1), M.sum(1)), "D_lower": np.nanquantile(boots, alpha, axis=0),
                             "D_upper": np.nanquantile(boots, 1 - alpha, axis=0)}, index=tab.index)


def _as_set(x) -> set:
    if isinstance(x, str):
        return {x}
    if isinstance(x, Iterable):
        return {str(v) for v in x}
    raise TypeError("group labels must be a string or a list of strings")
