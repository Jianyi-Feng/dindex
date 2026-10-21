from pathlib import Path
import numpy as np
import pandas as pd
import pytest
import dindex
from dindex import d_index, d_from_table, d_ci_counts, d_ci_bootstrap, log_ratio, d_from_log_ratio
from dindex.cli import main

DATA = Path(__file__).resolve().parents[1] / "examples" / "data"


def test_formula_and_bounds():
    assert d_index(3, 1) == pytest.approx(0.5)
    assert d_index(0, 5) == -1 and d_index(5, 0) == 1 and d_index(2, 2) == 0
    assert np.isnan(d_index(0, 0))


def test_scale_invariance_and_log_ratio_identity():
    p, m = np.array([0.3, 12.0, 7e-4]), np.array([0.1, 3.0, 2e-3])
    assert np.allclose(d_index(p, m), d_index(1000 * p, 1000 * m))
    assert np.allclose(d_index(p, m), d_from_log_ratio(log_ratio(p, m)))


def test_negative_rejected():
    with pytest.raises(ValueError):
        d_index(-1, 2)
    with pytest.raises(ValueError):
        log_ratio(-1, 2)
    with pytest.raises(ValueError):
        log_ratio(1, 2, pseudocount=0)


def test_table_groups_qc_and_overlap():
    ab = pd.DataFrame({"f1": [10, 0, 0], "f2": [5, 20, 0], "f3": [1, 1, 50], "f4": [4, 0, 0]}, index=["s1", "s2", "s3"])
    tax = {"f1": "A", "f2": "B", "f3": "C", "f4": "A"}
    r = d_from_table(ab, tax, plus="A", minus="B", min_share=0.5)
    assert r.loc["s1", "D"] == pytest.approx((14 - 5) / 19)
    assert r.loc["s1", "share"] == pytest.approx(19 / 20)
    assert r.loc["s3", "qc"] == "undefined" and np.isnan(r.loc["s3", "D"])
    assert r.loc["s2", "qc"] == "ok"
    r2 = d_from_table(ab.T, tax, plus=["A"], minus=["B", "C"], features_as="rows")
    assert r2.loc["s3", "D"] == -1
    with pytest.raises(ValueError):
        d_from_table(ab, tax, plus=["A", "B"], minus="B")
    with pytest.raises(ValueError):
        d_from_table(ab, tax, plus="A", minus="Z")


def test_count_interval_contains_estimate_and_narrows():
    est, lo, hi = d_ci_counts([30, 3000], [10, 1000])
    assert np.all(lo <= est) and np.all(est <= hi)
    assert (hi[1] - lo[1]) < (hi[0] - lo[0])
    with pytest.raises(ValueError):
        d_ci_counts([-1], [2])
    with pytest.raises(ValueError):
        d_ci_counts([1], [2], level=1)


def test_bootstrap_interval():
    rng = np.random.default_rng(1)
    ab = pd.DataFrame(rng.gamma(1, 1, (4, 30)), columns=[f"f{i}" for i in range(30)])
    tax = {f"f{i}": ("A" if i < 15 else "B") for i in range(30)}
    ci = d_ci_bootstrap(ab, tax, "A", "B", n_boot=300)
    assert ((ci.D_lower <= ci.D) & (ci.D <= ci.D_upper)).all()
    with pytest.raises(ValueError):
        d_ci_bootstrap(ab, tax, "A", "missing", n_boot=10)


def test_reproduces_published_spot_and_tara_values():
    spot = pd.read_csv(DATA / "SPOT_sample_Dindex.csv")
    assert np.allclose(d_index(spot.Imitervirales_abundance, spot.Algavirales_abundance), spot.D_index, atol=1e-9)
    tara = pd.read_csv(DATA / "Tara_observed_Dindex.csv")
    ok = tara.D_index_qc == "ok"
    assert np.allclose(d_index(tara.I_relative[ok], tara.A_relative[ok]), tara.D_index[ok], atol=1e-9)


def test_calibration_recovers_known_slope_and_zero_point():
    cal = pytest.importorskip("dindex.calibrate")
    rng = np.random.default_rng(3)
    t = rng.uniform(5, 30, 400)
    df = pd.DataFrame({"T": t, "D": np.clip(0.08 * (t - 17) + rng.normal(0, 0.05, 400), -1, 1)})
    s = cal.sensitivity(df, "T")
    assert s["ci_low"] < 0.08 < s["ci_high"]
    z = cal.zero_point(df, "T", n_boot=200)
    assert z["ci_low"] < 17 < z["ci_high"]


def test_seasonal_baseline_residuals():
    cal = pytest.importorskip("dindex.calibrate")
    dates = pd.date_range("2010-01-15", periods=72, freq="MS")
    d = 0.5 * np.sin(2 * np.pi * dates.dayofyear.to_numpy() / 365.25)
    d[60:] += 0.4
    res = cal.seasonal_baseline(pd.DataFrame({"date": dates, "D": d}), "date", baseline_end="2014-12-31")
    assert res.D_residual.iloc[:60].abs().max() < 0.05 and res.D_residual.iloc[60:].mean() > 0.3


def test_cli_compute(tmp_path):
    ab = pd.DataFrame({"f1": [10, 1], "f2": [5, 20], "f3": [1, 1]}, index=pd.Index(["s1", "s2"], name="sample"))
    tax = pd.DataFrame({"feature": ["f1", "f2", "f3"], "order": ["Imitervirales", "Algavirales", "Other"]})
    ab.to_csv(tmp_path / "ab.tsv", sep="\t"); tax.to_csv(tmp_path / "tax.tsv", sep="\t", index=False)
    out = tmp_path / "d.tsv"
    main(["compute", "--abundance", str(tmp_path / "ab.tsv"), "--taxonomy", str(tmp_path / "tax.tsv"), "--rank", "order",
          "--plus", "Imitervirales", "--minus", "Algavirales", "--counts", "-o", str(out)])
    r = pd.read_csv(out, sep="\t", index_col=0)
    assert r.loc["s1", "D"] == pytest.approx(5 / 15) and {"D_lower", "D_upper", "share", "qc"} <= set(r.columns)


def test_screen_pairs_finds_opposing_groups():
    cal = pytest.importorskip("dindex.calibrate")
    rng = np.random.default_rng(5)
    t = pd.Series(rng.uniform(5, 30, 60))
    g = pd.DataFrame({"warm": np.exp(0.15 * t) * rng.gamma(5, 1, 60), "cold": np.exp(-0.15 * t) * 400 * rng.gamma(5, 1, 60),
                      "neutral": 50 * rng.gamma(5, 1, 60)})
    top = cal.screen_pairs(g, t).iloc[0]
    assert (top.plus, top.minus) == ("warm", "cold")
    empty = cal.screen_pairs(g[["neutral"]], t)
    assert empty.empty and list(empty.columns) == ["plus", "minus", "correlation", "n", "combined_mean_share"]


def test_calibration_input_checks():
    cal = pytest.importorskip("dindex.calibrate")
    df = pd.DataFrame({"T": [1, 2, 3], "D": [-0.2, 0, 0.2], "station": ["a", "b", "b"]})
    with pytest.raises(ValueError):
        cal.sensitivity(df, "T", cluster="station", hac_lags=1)
    with pytest.raises(ValueError):
        cal.zero_point(df.iloc[:2], "T")
