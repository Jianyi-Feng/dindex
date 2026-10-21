"""Calibration and change-detection helpers for a D-index.

D is best read as change relative to a local baseline: its sensitivity to an environmental driver can transfer between systems
while its zero point often does not. These helpers estimate both and quantify departures from a seasonal baseline.
Requires statsmodels.
"""
from __future__ import annotations

from typing import Sequence

import numpy as np
import pandas as pd

__all__ = ["sensitivity", "zero_point", "seasonal_baseline", "screen_pairs"]


def _smf():
    try:
        import statsmodels.formula.api as smf
    except ImportError as e:  # pragma: no cover
        raise ImportError("dindex.calibrate requires statsmodels (pip install statsmodels)") from e
    return smf


def sensitivity(data: pd.DataFrame, driver: str, d_col: str = "D", covariates: Sequence[str] = (), cluster: str | None = None,
                hac_lags: int | None = None) -> dict:
    """Slope of D on an environmental driver (e.g. temperature), with optional categorical covariates.

    Errors: cluster-robust if `cluster` is given (e.g. station), Newey-West if `hac_lags` is given (time-ordered data),
    otherwise HC3. Returns slope, 95% CI, P and n.
    """
    if cluster is not None and hac_lags is not None:
        raise ValueError("choose either cluster-robust or HAC errors, not both")
    if hac_lags is not None and hac_lags < 0:
        raise ValueError("hac_lags must be non-negative")
    smf = _smf()
    cols = [d_col, driver, *covariates] + ([cluster] if cluster else [])
    df = data[cols].dropna()
    if len(df) < 3 or df[driver].nunique() < 2:
        raise ValueError("at least three complete observations and two driver values are required")
    if cluster is not None and df[cluster].nunique() < 2:
        raise ValueError("cluster-robust errors require at least two clusters")
    rhs = " + ".join([f"Q('{driver}')"] + [f"C(Q('{c}'))" for c in covariates])
    model = smf.ols(f"Q('{d_col}') ~ {rhs}", data=df)
    if cluster:
        fit = model.fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(df[cluster])[0]})
    elif hac_lags is not None:
        fit = model.fit(cov_type="HAC", cov_kwds={"maxlags": hac_lags})
    else:
        fit = model.fit(cov_type="HC3")
    term = f"Q('{driver}')"
    lo, hi = fit.conf_int().loc[term]
    return {"slope": float(fit.params[term]), "ci_low": float(lo), "ci_high": float(hi), "P": float(fit.pvalues[term]), "n": int(fit.nobs)}


def zero_point(data: pd.DataFrame, driver: str, d_col: str = "D", n_boot: int = 2000, group: str | None = None, seed: int = 0) -> dict:
    """Driver value at which a linear fit of D crosses zero (equal abundance of the two groups), with a bootstrap 95% CI.

    If `group` is given (e.g. station), whole groups are resampled.
    """
    df = data[[d_col, driver] + ([group] if group else [])].dropna()
    if not isinstance(n_boot, (int, np.integer)) or n_boot < 2:
        raise ValueError("n_boot must be an integer of at least 2")
    if len(df) < 3 or df[driver].nunique() < 2:
        raise ValueError("at least three complete observations and two driver values are required")
    if group is not None and df[group].nunique() < 2:
        raise ValueError("group bootstrap requires at least two groups")

    def x0(d):
        b1, b0 = np.polyfit(d[driver].to_numpy(float), d[d_col].to_numpy(float), 1)
        return -b0 / b1 if b1 != 0 else np.nan

    est = x0(df)
    rng = np.random.default_rng(seed)
    boots = []
    units = df[group].unique() if group else None
    for _ in range(n_boot):
        if group:
            pick = rng.choice(units, len(units))
            b = pd.concat([df[df[group] == u] for u in pick])
        else:
            b = df.sample(len(df), replace=True, random_state=int(rng.integers(1e9)))
        boots.append(x0(b))
    return {"zero_point": float(est), "ci_low": float(np.nanpercentile(boots, 2.5)), "ci_high": float(np.nanpercentile(boots, 97.5)), "n": len(df)}


def seasonal_baseline(data: pd.DataFrame, date_col: str, d_col: str = "D", baseline_end: str | None = None, harmonics: int = 1) -> pd.DataFrame:
    """Fit a seasonal harmonic model of D to a baseline period and return expected D and residuals for all samples.

    D ~ sum_k [sin(2 pi k doy/365.25) + cos(2 pi k doy/365.25)], fitted only to samples dated on or before `baseline_end`
    (all samples if None). Residuals (observed - expected) quantify departures such as marine heatwaves; analyse them at the
    level of years or events, because monthly residuals are autocorrelated.
    """
    if not isinstance(harmonics, (int, np.integer)) or harmonics < 1:
        raise ValueError("harmonics must be a positive integer")
    smf = _smf()
    df = data.copy()
    t = pd.to_datetime(df[date_col])
    doy = t.dt.dayofyear.to_numpy(float)
    terms = []
    for k in range(1, harmonics + 1):
        df[f"_s{k}"] = np.sin(2 * np.pi * k * doy / 365.25)
        df[f"_c{k}"] = np.cos(2 * np.pi * k * doy / 365.25)
        terms += [f"_s{k}", f"_c{k}"]
    base = df[(t <= pd.Timestamp(baseline_end)) if baseline_end else np.ones(len(df), bool)].dropna(subset=[d_col])
    if len(base) <= 2 * harmonics + 1:
        raise ValueError("baseline period has too few complete observations for the requested harmonics")
    fit = smf.ols(f"Q('{d_col}') ~ " + " + ".join(terms), data=base).fit()
    df["D_expected"] = fit.predict(df)
    df["D_residual"] = df[d_col] - df["D_expected"]
    df.attrs["baseline_r2"] = float(fit.rsquared)
    df.attrs["baseline_n"] = int(fit.nobs)
    return df.drop(columns=terms)


def screen_pairs(group_abundance: pd.DataFrame, driver: pd.Series, min_mean_share: float = 0.05, method: str = "spearman") -> pd.DataFrame:
    """Rank candidate group pairs by how strongly their D-index covaries with a driver.

    `group_abundance` is samples x groups (e.g. from summing a feature table by order). Only groups whose mean share of the
    table is at least `min_mean_share` are considered. Screening many pairs inflates false positives: choose the pair on a
    derivation dataset and evaluate it on independent data before interpreting it.
    """
    if not 0 <= min_mean_share <= 1:
        raise ValueError("min_mean_share must be between 0 and 1")
    vals = group_abundance.apply(pd.to_numeric, errors="raise")
    if (vals < 0).any().any():
        raise ValueError("group abundances must be non-negative")
    tot = vals.sum(axis=1).replace(0, np.nan)
    share = vals.div(tot, axis=0).mean()
    groups = share.index[share >= min_mean_share]
    rows = []
    for i, g1 in enumerate(groups):
        for g2 in groups[i + 1:]:
            d = pd.Series(np.where((vals[g1] + vals[g2]) > 0,
                                   (vals[g1] - vals[g2]) / (vals[g1] + vals[g2]).replace(0, np.nan),
                                   np.nan), index=vals.index)
            ok = d.notna() & driver.reindex(d.index).notna()
            if ok.sum() < 5:
                continue
            r = d[ok].corr(driver.reindex(d.index)[ok], method=method)
            plus, minus = (g1, g2) if r >= 0 else (g2, g1)
            rows.append({"plus": plus, "minus": minus, "correlation": abs(r), "n": int(ok.sum()),
                         "combined_mean_share": float(share[g1] + share[g2])})
    columns = ["plus", "minus", "correlation", "n", "combined_mean_share"]
    if not rows:
        return pd.DataFrame(columns=columns)
    return pd.DataFrame(rows, columns=columns).sort_values("correlation", ascending=False).reset_index(drop=True)
