"""Supplementary Table 12 workbook: final  CMIP6 results and the archived-model audit."""
from pathlib import Path
import pandas as pd
H = Path(__file__).resolve().parent; S = H.parents[1] / "processed_tables" / "supplementary_tables"
with pd.ExcelWriter(S / "Supplementary_Table_10_CMIP6_fronts.xlsx") as xw:
    pd.DataFrame({"sheet": ["A_rates_by_design", "B_front_trajectories", "C_agreement_by_band", "D_area_of_applicability", "E_predictor_units", "F_vs", "G_source_table_audit"],
                  "contents": ["Per-ESM poleward front rate (km per decade), R2, P and 2000s/2090s front latitude for each projection design",
                               "Front latitude by ESM, hemisphere and decade for each projection design (including calibration-domain-only fronts)",
                               "Area-weighted share of 2-degree cells where all four ESMs agree on the sign of dD (2090s minus 2000s), by absolute-latitude band and hemisphere",
                               "Area-weighted share of ocean within the Tara area of applicability, by predictor specification, ESM and decade",
                               "Co-located input comparison; the CO3 row is diagnostic only because the Tara field is a mislabelled pigment and was excluded",
                               "Cellwise comparison of the final seven-predictor  and archived eight-field  projections",
                               "Registry-level proof that the Tara field labelled carbonate is identical to HPLC 19'-butanoyloxyfucoxanthin"]}).to_excel(xw, sheet_name="README", index=False)
    pd.read_csv(H / "cmip6_front_rates_by_projection.csv").to_excel(xw, sheet_name="A_rates_by_design", index=False)
    pd.read_csv(H / "cmip6_front_trajectories_by_projection.csv").to_excel(xw, sheet_name="B_front_trajectories", index=False)
    pd.read_csv(H / "cmip6_agreement_by_band.csv").to_excel(xw, sheet_name="C_agreement_by_band", index=False)
    pd.read_csv(H / "cmip6_area_of_applicability.csv").to_excel(xw, sheet_name="D_area_of_applicability", index=False)
    pd.read_csv(H / "cmip6_vs_tara_predictor_units.csv").to_excel(xw, sheet_name="E_predictor_units", index=False)
    pd.read_csv(H / "cmip6_field_sensitivity_cellwise.csv").to_excel(xw, sheet_name="F_vs", index=False)
    pd.read_csv(H / "tara_carbonate_table_audit.csv").to_excel(xw, sheet_name="G_source_table_audit", index=False)
print("ok")
