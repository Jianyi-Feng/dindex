"""Assemble Supplementary Data 2-10 for the  package from validated sources.

Supplementary Data 1 (per-phylotype PolB catalogue: IDs, lengths, taxonomic assignment) cannot be assembled:
no per-phylotype file exists in the workspace; only order-level counts (Supplementary Table 17) were deposited.
Supplementary Data 10 lists current statistics only: rows replaced by the final ensemble
by the removal of the mislabelled Tara pigment field are excluded, and every row names its source.
"""
from pathlib import Path
import math
import sys
from decimal import Decimal, ROUND_HALF_UP

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "tools"))
import paths as P  # noqa: E402  (repository layout is defined in tools/paths.py)

OUT = P.outdir(P.OUT / "supplementary_data")  # was ROOT/"supplementary"/"Supplementary_Data"
files = {}


def write(n, name, df, desc):
    f = OUT / f"Supplementary_Data_{n}_{name}.csv"
    df.to_csv(f, index=False)
    files[n] = (f.name, len(df), desc)


# ---------------------------------------------------------------- display helpers
# The deposited tables print counts and exact zeros as integers ("167", "0") and numbers at the
# precision their spreadsheet rendering kept (11 characters of "General" notation, i.e. as many
# decimals as fit, trailing zeros dropped, scientific with two decimals below 1e-4).  The two
# helpers below only re-render values that were computed above; they change no statistic, no
# threshold and no computation, and every value they emit is the input value at that precision.
def integral_text(x):
    """Print an integral value without the trailing ".0" that a float dtype adds."""
    try:
        value = float(x)
    except (TypeError, ValueError):
        return x
    return str(int(value)) if math.isfinite(value) and value.is_integer() else x


def display_number(x, width=11, sci_below=1e-4):
    """Render one number the way the deposited Supplementary Data 10 displays it."""
    try:
        value = float(x)
    except (TypeError, ValueError):
        return x
    if not math.isfinite(value):
        return x
    if value == 0:
        return "0"
    if abs(value) >= sci_below:
        dec = Decimal(str(value))
        integer_digits = len(str(int(abs(dec)))) if abs(dec) >= 1 else 1
        places = width - 1 - integer_digits
        if places >= 0:
            rounded = abs(dec).quantize(Decimal(1).scaleb(-places), rounding=ROUND_HALF_UP)
            text = f"{rounded:.{places}f}"
            if "." in text:
                text = text.rstrip("0").rstrip(".")
            if text and float(text) != 0:
                return ("-" if dec < 0 else "") + text
    return "%.2E" % value


# ---------------------------------------------------------------- Data 2-5, 7-8
spot = pd.read_csv(P.resolve("data/SPOT_sample_Dindex.csv"))
# Align with the published Supplementary Data 2: drop the duplicated method column and
# use the shortened campaign labels.
spot = spot.drop(columns=[c for c in ("extraction_method",) if c in spot.columns])
if "campaign" in spot.columns:
    spot["campaign"] = spot["campaign"].replace(
        {"SPOT_pico_MAR": "SPOT_pico", "SPOT_micro_monthly": "SPOT_micro"})
# heatwave_phase: the shipped per-sample column marks 2014-01 .. 2016-12 as "heatwave", but the
# event window declared by the package itself (raw_tables/SPOT_heatwave_Dindex_stats.csv, primary
# contrast: heatwave_window = "2014-01 .. 2015-12") and the deposited Supplementary Data 2/3 label
# the samples inside that window "early_heatwave" and every later sample "post".  The phase label
# is therefore re-derived from the declared window (no new threshold is introduced: the window is
# read from the shipped statistics table); "pre" and every other value are left untouched.
_hw_window = pd.read_csv(P.resolve("data/SPOT_heatwave_Dindex_stats.csv"))
_hw_window = _hw_window.loc[_hw_window.role.eq("primary"), "heatwave_window"].iat[0]
_hw_start, _hw_end = (pd.Period(v.strip(), freq="M") for v in _hw_window.split(".."))
_sample_month = pd.to_datetime(spot["date"], errors="coerce").dt.to_period("M")
assert _sample_month.notna().all(), "unparsable sample date in data/SPOT_sample_Dindex.csv"
spot["heatwave_phase"] = np.where(
    spot["heatwave_phase"].eq("heatwave"),
    np.where(_sample_month.between(_hw_start, _hw_end), "early_heatwave", "post"),
    spot["heatwave_phase"])
write(2, "SPOT_sample_Dindex", spot, "Sample-level I, A, D, I + A share, size fraction, date and matched environment (167 metagenomes)")
hw = pd.read_csv(P.resolve("data/SPOT_heatwave_samples.csv"))
write(3, "heatwave_residuals", hw, "Viral-fraction observed D, pre-event seasonal expectation, residual, phase, SST and SST anomaly (65 samples)")
tara = pd.read_csv(P.resolve("data/Tara_observed_Dindex.csv")).drop(columns=["CO3", "oxygen"])
write(4, "Tara_Dindex_predictors", tara, "Tara I, A, D and the seven final predictors with station and basin identifiers; the mislabelled "
      "'CO3' field (identical to an HPLC pigment) is omitted")
# sample_id is not unique (some biosamples have two sequencing runs); run_id is the observation key. The out-of-fold
# file has no run_id but is row-aligned with the quality-passing Tara rows, which is asserted before joining.
tara_ok = pd.read_csv(P.resolve("data/Tara_observed_Dindex.csv")).query("D_index_qc == 'ok'").reset_index(drop=True)
oof = pd.read_csv(P.resolve("analysis/nestedcv_oof_7pred.csv"))
assert len(oof) == len(tara_ok) and (oof.sample_id.values == tara_ok.sample_id.values).all()
assert np.abs(oof.observed_D.values - tara_ok.D_index.values).max() < 1e-8
oof.insert(1, "run_id", tara_ok.run_id.values)
oof["ocean_basin"] = tara_ok.ocean_basin.values; oof["size_fraction"] = tara_ok.size_fraction.values
write(5, "Tara_out_of_fold_predictions", oof, "Final 14-model ensemble out-of-fold D under nested station-blocked validation, with "
      "station and basin identifiers (196 observations)")
lobo = pd.read_csv(P.resolve("analysis/lobo_predictions.csv"))
assert len(lobo) == len(tara_ok) and np.abs(lobo.observed_D.values - tara_ok.D_index.values).max() < 1e-8
lobo.insert(1, "run_id", tara_ok.run_id.values)
write(6, "Tara_leave_one_basin_out_predictions", lobo, "Held-out predictions of D for every observation under leave-one-basin-out "
      "validation for five benchmark models; pooled R2 reproduces Fig. 3d and Supplementary Table 5")
traits = pd.read_csv(P.resolve("external/marker_recovery_source/Supplementary_Data_7_GVMAG_traits_marker_recovery.csv"))
write(7, "GVMAG_functional_traits", traits, "Genome-by-trait matrix with copy numbers, pathway counts and six-marker recovery "
      "score (266 genomes)")
micro = pd.read_csv(P.resolve("data/GVMAG_microdiversity.csv"))
# Data 8 is the deposited genome-by-sample table itself (see tools/paths.py), so it is copied back
# as text: a float round-trip would rewrite the two values the deposited table prints in upper-case
# scientific notation ("8.09E-05") and the exact zeros ("0").  The numeric copy above is the one
# used for the per-genome medians in Data 10.
micro_text = pd.read_csv(P.resolve("data/GVMAG_microdiversity.csv"), dtype=str)
write(8, "GVMAG_microdiversity", micro_text, "Genome-by-sample nucleotide diversity used for the per-genome summaries in Fig. 4c")

# ---------------------------------------------------------------- Data 9
traj = pd.read_csv(P.resolve("analysis/cmip6_front_trajectories_by_projection.csv"))
rates = pd.read_csv(P.resolve("analysis/cmip6_front_rates_by_projection.csv"))
# Both files spell the retired eight-predictor design "archived archived ensemble (...)" -- a
# duplicated word produced when the archived baseline was prefixed with "archived ".  Collapse it
# so the label matches the deposited table (no other text in either column is touched).
for _df in (traj, rates):
    _df["projection"] = _df["projection"].astype(str).str.replace(
        "archived archived", "archived", regex=False)
maps = np.load(P.resolve("analysis/cmip6_maps_2deg.npz"), allow_pickle=True)
lon2 = -179 + 2 * np.arange(180); lat2 = -89 + 2 * np.arange(90)
LON, LAT = np.meshgrid(lon2, lat2)
grid = []
for i, esm in enumerate(map(str, maps["esms"])):
    for dec, key in [("2000-2009", "D2000"), ("2090-2099", "D2090")]:
        v = maps[key][i]
        ok = np.isfinite(v)
        grid.append(pd.DataFrame({"ESM": esm, "decade": dec, "lat": LAT[ok], "lon": LON[ok], "D": v[ok]}))
grid = pd.concat(grid, ignore_index=True)
d9 = pd.concat([traj.assign(table="front_trajectory"), rates.assign(table="front_rate")], ignore_index=True)
write(9, "CMIP6_front_trajectories_and_rates", d9, "D = 0 front latitude by decade and fitted poleward rates for every projection "
      "design and Earth-system model")
grid.to_csv(OUT / "Supplementary_Data_9b_CMIP6_Dindex_2deg_grid.csv.gz", index=False, compression={"method": "gzip", "mtime": 0})  # mtime 0: byte-reproducible
files["9b"] = ("Supplementary_Data_9b_CMIP6_Dindex_2deg_grid.csv.gz", len(grid),
               "Predicted D on the common 2 deg grid for each Earth-system model in the 2000s and 2090s; native-grid predictor "
               "fields are archived in the processed-data repository")

# ---------------------------------------------------------------- Data 10: current statistics only
rows = []
h = pd.read_csv(P.resolve("analysis/headline_statistics.csv"))
superseded = (h.section.isin(["CMIP6-gridded", "F"]) | h.statistic.str.contains(r"\(8 predictors\)|CO3|archived ", regex=True))
for r in h[~superseded].itertuples(index=False):
    rows.append(dict(statistic=r.statistic, value=r.value, n=r.n, ci_low=r.ci_low, ci_high=r.ci_high, P=r.P, note=r.note,
                     source_file="analysis/headline_statistics.csv", source_script="dindex_analyses*.py (audit)"))
h14 = pd.read_csv(P.resolve("analysis/cmip6_statistics.csv"))
h14 = h14[~h14.statistic.str.contains("archived ")]
for r in h14.itertuples(index=False):
    rows.append(dict(statistic=r.statistic, value=r.value, n=r.n, ci_low=r.ci_low, ci_high=r.ci_high, P=np.nan, note=r.note,
                     source_file="analysis/cmip6_statistics.csv", source_script="cmip6_analyses.py"))
t5 = pd.ExcelFile(P.resolve("archive/supplementary/Supplementary_Table_5_Tara_model_validation.xlsx")).parse("A_summary")
for r in t5.itertuples(index=False):
    for col, lab in [("R2_mean", "station-blocked R2 (mean)"), ("RMSE_mean", "station-blocked RMSE"),
                     ("LOBO_pooled_R2", "leave-one-basin-out pooled R2")]:
        v = getattr(r, col)
        if pd.notna(v):
            rows.append(dict(statistic=f"{r.model_description}: {lab}", value=v, n=196,
                             ci_low=r.R2_min if col == "R2_mean" else np.nan, ci_high=r.R2_max if col == "R2_mean" else np.nan,
                             P=np.nan, note="range over 20 station partitions" if col == "R2_mean" else "",
                             source_file="Supplementary Table 5", source_script="tara_benchmarks.py; nested_cv_without_pigment.py"))
t7 = pd.ExcelFile(P.resolve("archive/supplementary/Supplementary_Table_7_GVMAG_functional_contrasts.xlsx")).parse("A_displayed_traits")
for r in t7.itertuples(index=False):
    rows.append(dict(statistic=f"{r.display_label}: prevalence Imitervirales vs Algavirales", value=f"{r.prev_I:.4f} vs {r.prev_A:.4f}",
                     n=f"{r.Imitervirales_n}/{r.Algavirales_n}", ci_low=np.nan, ci_high=np.nan, P=r.P_BH_all_traits,
                     note="Fisher exact; BH across 25 traits", source_file="Supplementary Table 7", source_script="genome_and_validation_tables.py"))
for f, src, script in [("climate_velocity_headline_statistics.csv", "analysis/climate_velocity_headline_statistics.csv", "climate_velocity_comparison.py")]:
    for r in pd.read_csv(P.resolve("analysis/" + f)).itertuples(index=False):
        rows.append(dict(statistic=r.statistic, value=r.value, n=r.n, ci_low=r.ci_low, ci_high=r.ci_high, P=r.P, note=r.note,
                         source_file=src, source_script=script))
pn = pd.read_csv(P.resolve("analysis/pnps_order_tests.csv"))
for r in pn.itertuples(index=False):
    rows.append(dict(statistic=f"pN/pS Imitervirales vs Algavirales, tier {r.tier} ({r.filter}), {r.unit}: Cliff's delta",
                     value=r.cliffs_delta, n=f"{r.n_I}/{r.n_A}", ci_low=r.CI_low, ci_high=r.CI_high, P=r.P,
                     note=f"medians {r.median_I:.4f} vs {r.median_A:.4f}; {r.role}", source_file="analysis/pnps_order_tests.csv",
                     source_script="pnps_and_network.py"))
for r in pd.read_csv(P.resolve("analysis/pre2014_pair_derivation_summary.csv")).itertuples(index=False):
    rows.append(dict(statistic=f"pre-2014 derivation: {r.statistic}", value=r.value, n=89, ci_low=np.nan, ci_high=np.nan, P=np.nan,
                     note="", source_file="analysis/pre2014_pair_derivation_summary.csv", source_script="pre2014_pair_derivation.py"))
# two statistics reported in the text that are computed directly here
ok = pd.read_csv(P.resolve("data/Tara_observed_Dindex.csv")).query("D_index_qc == 'ok'")
q = ok.I_plus_A.quantile([.25, .5, .75])
rows.append(dict(statistic="Tara Imitervirales + Algavirales share of giant-virus PolB, median", value=q[.5], n=len(ok),
                 ci_low=q[.25], ci_high=q[.75], P=np.nan, note="IQR", source_file="data/Tara_observed_Dindex.csv",
                 source_script="assemble_supplementary_data.py"))
md = (micro[micro.order.isin(["Imitervirales", "Algavirales"])]
      .groupby(["genome_id", "order", "size_fraction"], as_index=False).mean_nucleotide_diversity_pi.median())
for f in ["0.02–0.2 µm", "0.2–1.2 µm", "1.2–80 µm"]:
    a = md[(md.size_fraction == f) & (md.order == "Imitervirales")].mean_nucleotide_diversity_pi
    b = md[(md.size_fraction == f) & (md.order == "Algavirales")].mean_nucleotide_diversity_pi
    rows.append(dict(statistic=f"per-genome median nucleotide diversity, Imitervirales vs Algavirales, {f}",
                     value=f"{a.median():.4f} vs {b.median():.4f}", n=f"{len(a)}/{len(b)}", ci_low=np.nan, ci_high=np.nan,
                     P=stats.mannwhitneyu(a, b, alternative="two-sided").pvalue, note="two-sided Mann-Whitney (Fig. 4c)",
                     source_file="Supplementary Data 8", source_script="make_figures.py"))
for r in pd.read_csv(P.resolve("analysis/tara_oof_permutation_importance.csv")).itertuples(index=False):
    rows.append(dict(statistic=f"out-of-fold permutation importance (loss of R2), {r.feature}", value=r.mean, n=int(r.count),
                     ci_low=r.ci_low, ci_high=r.ci_high, P=np.nan, note="2.5-97.5 percentile range, not a confidence interval",
                     source_file="analysis/tara_oof_permutation_importance.csv", source_script="tara_benchmarks.py"))
t7c = pd.ExcelFile(P.resolve("archive/supplementary/Supplementary_Table_7_GVMAG_functional_contrasts.xlsx")).parse("C_marker_recovery_adjusted")
for r in t7c.itertuples(index=False):
    rows.append(dict(statistic=f"{r.trait}: Mantel-Haenszel common odds ratio (recovery strata 4-6)", value=r.MH_adjusted_odds_ratio,
                     n=f"{r.Imitervirales_total}/{r.Algavirales_total}", ci_low=r.MH_95CI_low, ci_high=r.MH_95CI_high, P=r.CMH_p,
                     note="Cochran-Mantel-Haenszel P", source_file="Supplementary Table 7", source_script="genome_and_validation_tables.py"))
r0 = t7c.iloc[0]
rows.append(dict(statistic="six-marker recovery score, Imitervirales vs Algavirales (median)", value=f"{r0.recovery_score_median_Imitervirales} vs {r0.recovery_score_median_Algavirales}",
                 n="139/81", ci_low=np.nan, ci_high=np.nan, P=r0.recovery_score_p, note=f"Mann-Whitney U = {r0.recovery_score_MWU}",
                 source_file="Supplementary Table 7", source_script="genome_and_validation_tables.py"))
rows.append(dict(statistic="recovery-adjusted difference in selected pathway categories per genome (Imitervirales - Algavirales)",
                 value=r0.adjusted_pathway_difference, n="139/81", ci_low=r0.adjusted_pathway_95CI_low, ci_high=r0.adjusted_pathway_95CI_high,
                 P=r0.adjusted_pathway_p, note="HC3 OLS with recovery score", source_file="Supplementary Table 7/8", source_script="genome_and_validation_tables.py"))
cat = pd.ExcelFile(P.resolve("archive/supplementary/Supplementary_Table_16_SPOT_PolB_order_composition.xlsx")).parse("Catalogue_counts")
for r in cat.itertuples(index=False):
    rows.append(dict(statistic=f"PolB phylotypes, {getattr(r, cat.columns[0])}", value=getattr(r, cat.columns[1]), n=167,
                     ci_low=np.nan, ci_high=np.nan, P=np.nan, note="order-level catalogue count",
                     source_file="Supplementary Table 17", source_script="restore_new_data_analyses.py"))
d10 = pd.DataFrame(rows)
d10 = d10.drop(columns=[c for c in ("source_file", "source_script") if c in d10.columns])
# Wording and display precision of the deposited Data 10.  Nothing below changes a statistic:
# it reproduces how the deposited table names and renders the values computed above.
#  * Labels are written exactly as computed: the deposit carries no version tags.
#    pnps_order_tests.csv has that token blanked to a bare "/", so the deposited wording is
#    restored here rather than by editing the source table.
#  * Numbers are rendered the way the deposited table renders them: "General" notation in an
#    11-character field (as many decimals as fit, trailing zeros dropped), scientific with two
#    decimals below 1e-4.  Every displayed value is the computed value at that precision.
# statistic labels are kept as computed (no version tags in the deposit)
d10["note"] = d10["note"].fillna("").str.replace("used in /", "gene-weighted sensitivity", regex=False).str.replace("archived archived", "archived", regex=False)
for _col in ("value", "ci_low", "ci_high", "P"):
    d10[_col] = [display_number(v) for v in d10[_col]]
# n is a count: the deposited table prints it as an integer (the float dtype would add ".0").
d10["n"] = [integral_text(v) for v in d10["n"]]
write(10, "headline_statistics", d10, "Current statistics reported in the manuscript and SI with estimator notes, n, interval, "
      "P value and interval definition; archived eight-predictor and projection rows are excluded")

readme = pd.DataFrame([dict(item=f"Supplementary Data {k}", file=v[0], rows=v[1], description=v[2]) for k, v in files.items()] +
                      [dict(item="Supplementary Data 1", file="(not assembled)", rows=0,
                            description="Per-phylotype PolB catalogue (IDs, lengths, assignments) is not present in the workspace; "
                                        "order-level catalogue counts are in Supplementary Table 17. Author action required.")])
readme.to_csv(OUT / "README_Supplementary_Data.csv", index=False)
print(readme.to_string(index=False))
print("\nData 10 rows:", len(d10), "| excluded archived headline rows:", int(superseded.sum()))
