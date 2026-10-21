"""Worked example: reproduce the core D-index analyses of the SPOT/Tara giant-virus study with the dindex toolkit.

Run from the toolkit root:  python examples/spot_tara_example.py
"""
from pathlib import Path
import pandas as pd
import dindex
from dindex import calibrate

DATA = Path(__file__).resolve().parent / "data"

# 1. SPOT: summed PolB coverage of the two orders -> D (identical to the published values)
spot = pd.read_csv(DATA / "SPOT_sample_Dindex.csv", parse_dates=["date"])
spot["D"] = dindex.d_index(spot.Imitervirales_abundance, spot.Algavirales_abundance)
spot["share"] = spot.I_plus_A
print("max |D - published D| =", float((spot.D - spot.D_index).abs().max()))
print("median share of PolB coverage in the two orders:", round(spot.share.median(), 3))

# 2. Temperature sensitivity with size fraction as a covariate, then per fraction with Newey-West errors
print("SPOT sensitivity:", calibrate.sensitivity(spot, "SST", covariates=["size_fraction_ascii"]))
for frac, g in spot.sort_values("date").groupby("size_fraction_ascii"):
    s = calibrate.sensitivity(g, "SST", hac_lags=3)
    z = calibrate.zero_point(g, "SST", n_boot=500)
    print(f"  {frac}: slope {s['slope']:.3f} ({s['ci_low']:.3f}-{s['ci_high']:.3f}); D = 0 at {z['zero_point']:.1f} degC "
          f"({z['ci_low']:.1f}-{z['ci_high']:.1f})")

# 3. Seasonal baseline on pre-2014 viral-fraction samples and the 2014-2015 marine-heatwave displacement
viral = spot[spot.size_fraction_ascii == "0.02-0.2um"]
res = calibrate.seasonal_baseline(viral, "date", baseline_end="2013-12-31")
pre = res[res.date < "2014-01-01"].D_residual.mean(); hw = res[(res.date >= "2014-01-01") & (res.date <= "2015-12-31")].D_residual.mean()
print(f"heatwave residual shift (2014-2015 minus 2009-2013): {hw - pre:+.2f} D units; baseline R2 = {res.attrs['baseline_r2']:.2f}")

# 4. Tara: station-clustered sensitivity and the zero point with station bootstrap
tara = pd.read_csv(DATA / "Tara_observed_Dindex.csv")
tara = tara[tara.D_index_qc == "ok"].copy()
tara["D"] = dindex.d_index(tara.I_relative, tara.A_relative)
print("Tara sensitivity (station-clustered):", calibrate.sensitivity(tara, "temperature", covariates=["size_fraction"], cluster="station_id"))
