# -*- coding: utf-8 -*-
"""Central path constants for the dindex reproducibility package.

Why this module exists
----------------------
The scripts in ``analysis/`` and ``figures/`` were written inside the author's
working tree, which used the flat layout::

    data/            analysis/         manuscript/
    supplementary/   external/         archive/supplementary/

The published repository ships a different (auditable) layout::

    raw_tables/      processed_tables/{supplementary_data, supplementary_tables,
                                      figure_source_data, model_and_statistics,
                                      phylogenies, cmip6_gridded_fields,
                                      genomes_and_microdiversity}

This module is the single place where that translation lives, so that no
script has to hard-code a directory again.  **It only maps paths.**  No number,
statistic, threshold, column selection or plotting parameter is defined here.

Usage
-----
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
    import paths as P

    df = pd.read_csv(P.RAW / "SPOT_sample_Dindex.csv")          # new-style
    df = pd.read_csv(P.resolve("data/SPOT_sample_Dindex.csv"))  # old-style name

``P.resolve`` takes the *historical* path exactly as the original script wrote
it and returns the file that is actually shipped.  It raises ``FileNotFoundError``
with the search record when the file is absent, so a missing input is reported
at the point of use instead of as a bare IOError.

Outputs
-------
Everything the scripts write goes under :data:`OUT` (default ``<repo>/outputs``,
git-ignored).  Set the environment variable ``DINDEX_OUT`` to redirect the
outputs, e.g. for a clean-room run::

    set DINDEX_OUT=D:\...\e2e_run\out

The ``*.npz`` gridded CMIP6 arrays are **not** part of the git repository
(they exceed the GitHub size limit); download them from the figshare deposit
(DOI 10.6084/m9.figshare.34419060, folder ``07_CMIP6_gridded_fields/``) and put
them in ``processed_tables/cmip6_gridded_fields/``.  :func:`require_gridded`
reports exactly which files are missing.
"""
from __future__ import annotations

import os
from pathlib import Path

# --------------------------------------------------------------------- root
PKG = Path(__file__).resolve().parents[1]

# --------------------------------------------------------------- input trees
RAW = PKG / "raw_tables"
PROC = PKG / "processed_tables"

SUPP_DATA = PROC / "supplementary_data"
SUPP_TABLES = PROC / "supplementary_tables"
FIG_SOURCE = PROC / "figure_source_data"
MODEL_STATS = PROC / "model_and_statistics"
PHYLO = PROC / "phylogenies"
CMIP6 = PROC / "cmip6_gridded_fields"
GENOMES = PROC / "genomes_and_microdiversity"
INSTRAIN = GENOMES / "instrain_all_genes"
INSTRAIN_PATHWAY = GENOMES / "instrain_focus_on_special_pathway"

DOCS = PKG / "docs"
EXAMPLES = PKG / "examples"
R_CODE = PKG / "R"
SRC_CODE = PKG / "src"

# ------------------------------------------------------------------- outputs
#: Root for every generated file.  Override with the DINDEX_OUT environment
#: variable to keep the checkout clean (used by the clean-room end-to-end run).
OUT = Path(os.environ.get("DINDEX_OUT") or (PKG / "outputs")).resolve()
OUT_SOURCE_DATA = OUT / "source_data"
OUT_FIGURES = OUT / "figures"
OUT_TABLES = OUT / "tables"

#: The six gridded CMIP6 arrays that live in the figshare deposit.
GRIDDED_FIELDS = (
    "cmip6_grid.npz",
    "archived__cmip6_grid.npz",
    "cmip6_maps_2deg.npz",
    "archived__cmip6_maps_2deg.npz",
    "cmip6_aoa_2deg.npz",
    "archived__cmip6_aoa_2deg.npz",
)

# ---------------------------------------------------------------- legacy map
#: Historical path (as written by the original scripts) -> shipped file.
#:
#: Three kinds of entries:
#:   * plain moves      -- ``data/X.csv`` now lives in ``raw_tables/X.csv``;
#:   * table renames    -- the working-tree copies ``analysis/X.csv`` were
#:     published under ``processed_tables/model_and_statistics/`` and, when a
#:     later revision superseded them, with an ``inherited__`` prefix that marks
#:     them as carried over unchanged from the previous submission;
#:   * deposit moves    -- ``archive/supplementary/*.xlsx`` were published as
#:     ``processed_tables/supplementary_tables/*.xlsx`` (Supplementary Table
#:     numbering is unchanged).
LEGACY = {
    # ---- raw observations ------------------------------------------------
    "data/SPOT_sample_Dindex.csv": RAW / "SPOT_sample_Dindex.csv",
    "data/Tara_observed_Dindex.csv": RAW / "Tara_observed_Dindex.csv",
    "data/SPOT_18S_GV_network_nodes.csv": RAW / "SPOT_18S_GV_network_nodes.csv",
    "data/SPOT_18S_GV_network_edges.csv": RAW / "SPOT_18S_GV_network_edges.csv",
    "data/Dindex_model_validation.csv": RAW / "Dindex_model_validation.csv",
    "data/SPOT_heatwave_Dindex_stats.csv": RAW / "SPOT_heatwave_Dindex_stats.csv",
    "data/CMIP6_Dindex_latitude_year.csv": RAW / "CMIP6_Dindex_latitude_year.csv",
    "data/Dindex_front_trajectory.csv": RAW / "Dindex_front_trajectory.csv",
    # ---- derived tables that became Supplementary Data ------------------
    # The vectorised heatwave table (65 samples) was published as Supplementary
    # Data 3; it is the same table the script rebuilds from it.
    "data/SPOT_heatwave_samples.csv":
        SUPP_DATA / "Supplementary_Data_3_heatwave_residuals.csv",
    "data/GVMAG_microdiversity.csv":
        SUPP_DATA / "Supplementary_Data_8_GVMAG_microdiversity.csv",
    "data/GVMAG_functional_traits.csv":
        SUPP_DATA / "Supplementary_Data_7_GVMAG_functional_traits.csv",
    # ---- inStrain pN/pS gene-by-sample table -----------------------------
    "data/pnps_network/instrain_all_genes/inStrain_allgenes_pNpS_gene_by_sample_defined.csv":
        INSTRAIN / "inStrain_allgenes_pNpS_gene_by_sample_defined.csv",
    # ---- marker-recovery tables (external/marker_recovery_source/) ------
    "external/marker_recovery_source/Supplementary_Data_7_GVMAG_traits_marker_recovery.csv":
        GENOMES / "Supplementary_Data_7_GVMAG_traits_marker_recovery.csv",
    "external/marker_recovery_source/marker_recovery_adjusted_genome_comparisons.csv":
        GENOMES / "marker_recovery_adjusted_genome_comparisons.csv",
    # ---- analysis outputs (working-tree analysis/*.csv) -----------------
    "analysis/headline_statistics.csv": MODEL_STATS / "headline_statistics.csv",
    "analysis/cmip6_statistics.csv": MODEL_STATS / "cmip6_statistics.csv",
    "analysis/tara_cross_validation.csv": MODEL_STATS / "tara_cross_validation.csv",
    "analysis/tara_oof_permutation_importance.csv":
        MODEL_STATS / "tara_oof_permutation_importance.csv",
    "analysis/nestedcv_oof_7pred.csv": MODEL_STATS / "nestedcv_oof_7pred.csv",
    "analysis/ensemble_validation.csv": MODEL_STATS / "ensemble_validation.csv",
    "analysis/lobo_predictions.csv": MODEL_STATS / "lobo_predictions.csv",
    "analysis/pnps_order_tests.csv": MODEL_STATS / "pnps_order_tests.csv",
    "analysis/pre2014_pair_derivation_summary.csv":
        MODEL_STATS / "pre2014_pair_derivation_summary.csv",
    "analysis/cmip6_front_rates_by_projection.csv":
        MODEL_STATS / "cmip6_front_rates_by_projection.csv",
    # the front trajectories were published with the gridded fields
    "analysis/cmip6_front_trajectories_by_projection.csv":
        CMIP6 / "cmip6_front_trajectories_by_projection.csv",
    "analysis/cmip6_zonal_profiles.csv": CMIP6 / "cmip6_zonal_profiles.csv",
    "analysis/cmip6_gridded_zonal_profiles.csv": CMIP6 / "cmip6_gridded_zonal_profiles.csv",
    "analysis/cmip6_agreement_by_band.csv": MODEL_STATS / "cmip6_agreement_by_band.csv",
    "analysis/cmip6_vs_tara_predictor_units.csv":
        MODEL_STATS / "cmip6_vs_tara_predictor_units.csv",
    "analysis/cmip6_area_of_applicability.csv":
        MODEL_STATS / "cmip6_area_of_applicability.csv",
    "analysis/cmip6_front_vs_isotherm_by_model.csv":
        MODEL_STATS / "cmip6_front_vs_isotherm_by_model.csv",
    "analysis/cmip6_isotherm_displacement.csv":
        MODEL_STATS / "cmip6_isotherm_displacement.csv",
    "analysis/climate_velocity_headline_statistics.csv":
        MODEL_STATS / "climate_velocity_headline_statistics.csv",
    # ---- analysis outputs carried over unchanged (inherited__) ----------
    "analysis/spot_viral_heatwave_residuals.csv":
        MODEL_STATS / "inherited__spot_viral_heatwave_residuals.csv",
    "analysis/spot_transfer_from_tara.csv":
        MODEL_STATS / "inherited__spot_transfer_from_tara.csv",
    "analysis/tara_station_means.csv":
        MODEL_STATS / "inherited__tara_station_means.csv",
    "analysis/cmip6_front_per_esm.csv":
        MODEL_STATS / "inherited__cmip6_front_per_esm.csv",
    # ---- gridded CMIP6 arrays (figshare deposit) ------------------------
    "analysis/cmip6_maps_2deg.npz": CMIP6 / "cmip6_maps_2deg.npz",
    "analysis/archived__cmip6_maps_2deg.npz": CMIP6 / "archived__cmip6_maps_2deg.npz",
    "analysis/cmip6_aoa_2deg.npz": CMIP6 / "cmip6_aoa_2deg.npz",
    "analysis/archived__cmip6_aoa_2deg.npz": CMIP6 / "archived__cmip6_aoa_2deg.npz",
    "analysis/cmip6_grid.npz": CMIP6 / "cmip6_grid.npz",
    "analysis/archived__cmip6_grid.npz": CMIP6 / "archived__cmip6_grid.npz",
    "data/pnps_network/cmip6_maps_2deg.npz": CMIP6 / "cmip6_maps_2deg.npz",
    "data/pnps_network/cmip6_aoa_2deg.npz": CMIP6 / "cmip6_aoa_2deg.npz",
    "data/pnps_network/cmip6_grid.npz": CMIP6 / "cmip6_grid.npz",
    # ---- archive/supplementary -> processed_tables/supplementary_tables --
    "archive/supplementary/Supplementary_Table_1_SPOT_samples.xlsx":
        SUPP_TABLES / "Supplementary_Table_1_SPOT_samples.xlsx",
    "archive/supplementary/Supplementary_Table_2_Dindex_calibration.xlsx":
        SUPP_TABLES / "Supplementary_Table_2_Dindex_calibration.xlsx",
    "archive/supplementary/Supplementary_Table_3_heatwave_tests.xlsx":
        SUPP_TABLES / "Supplementary_Table_3_heatwave_tests.xlsx",
    "archive/supplementary/Supplementary_Table_4_Tara_samples.xlsx":
        SUPP_TABLES / "Supplementary_Table_4_Tara_samples.xlsx",
    "archive/supplementary/Supplementary_Table_5_Tara_model_validation.xlsx":
        SUPP_TABLES / "Supplementary_Table_5_Tara_model_validation.xlsx",
    "archive/supplementary/Supplementary_Table_6_permutation_importance.xlsx":
        SUPP_TABLES / "Supplementary_Table_6_permutation_importance.xlsx",
    "archive/supplementary/Supplementary_Table_7_GVMAG_functional_contrasts.xlsx":
        SUPP_TABLES / "Supplementary_Table_7_GVMAG_functional_contrasts.xlsx",
    "archive/supplementary/Supplementary_Table_8_GVMAG_functional_repertoire.xlsx":
        SUPP_TABLES / "Supplementary_Table_8_GVMAG_functional_repertoire.xlsx",
    "archive/supplementary/Supplementary_Table_9_GVMAG_microdiversity.xlsx":
        SUPP_TABLES / "Supplementary_Table_9_GVMAG_microdiversity.xlsx",
    "archive/supplementary/Supplementary_Table_10_CMIP6_fronts.xlsx":
        SUPP_TABLES / "Supplementary_Table_10_CMIP6_fronts.xlsx",
    "archive/supplementary/Supplementary_Table_11_software_databases_parameters.xlsx":
        SUPP_TABLES / "Supplementary_Table_11_software_databases_parameters.xlsx",
    "archive/supplementary/Supplementary_Table_12_gene_glossary.xlsx":
        SUPP_TABLES / "Supplementary_Table_12_gene_glossary.xlsx",
    "archive/supplementary/Supplementary_Table_13_isotherm_and_climate_velocity.xlsx":
        SUPP_TABLES / "Supplementary_Table_13_isotherm_and_climate_velocity.xlsx",
    "archive/supplementary/Supplementary_Table_14_pNpS.xlsx":
        SUPP_TABLES / "Supplementary_Table_14_pNpS.xlsx",
    "archive/supplementary/Supplementary_Table_15_association_network.xlsx":
        SUPP_TABLES / "Supplementary_Table_15_association_network.xlsx",
    "archive/supplementary/Supplementary_Table_16_SPOT_PolB_order_composition.xlsx":
        SUPP_TABLES / "Supplementary_Table_16_SPOT_PolB_order_composition.xlsx",
    "archive/supplementary/Supplementary_Table_17_phylogeny_metadata.xlsx":
        SUPP_TABLES / "Supplementary_Table_17_phylogeny_metadata.xlsx",
    "archive/Supplementary_Table_14_pNpS.xlsx":
        SUPP_TABLES / "Supplementary_Table_14_pNpS.xlsx",
}

#: Historical output directories -> where they are written now.
#: ``supplementary/Supplementary_Data`` etc. no longer exist in the repository;
#: generated tables go to OUT so that a run never overwrites delivered data.
LEGACY_OUT = {
    "supplementary/Supplementary_Data": OUT / "supplementary_data",
    "supplementary/figures": OUT_FIGURES,
    "supplementary": OUT_FIGURES,
    "archive/supplementary": OUT_TABLES,
    "out": OUT,
    "out/supplementary": OUT_FIGURES,
}


def resolve(name: str) -> Path:
    """Return the shipped file for a historical path, or raise with a hint."""
    key = str(name).replace("\\", "/")
    try:
        path = LEGACY[key]
    except KeyError:
        hint = "\n  ".join(sorted(LEGACY))
        raise KeyError(
            f"{name!r} is not in paths.LEGACY; add the mapping to tools/paths.py.\n"
            f"Known historical inputs:\n  {hint}"
        ) from None
    if not path.exists():
        raise FileNotFoundError(
            f"{name} -> {path} (missing).\n"
            f"The file is part of the repository; if you are running from a partial "
            f"checkout, restore {path.parent} from the deposit (figshare DOI "
            f"10.6084/m9.figshare.34419060)."
        )
    return path


def outdir(path) -> Path:
    """Create and return an output directory (idempotent)."""
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def missing_gridded() -> list:
    """Return the gridded CMIP6 arrays that are not present in the checkout."""
    return [n for n in GRIDDED_FIELDS if not (CMIP6 / n).exists()]


def require_gridded(*names: str) -> list:
    """Return the requested gridded arrays, raising if any is absent."""
    absent = [n for n in names if not (CMIP6 / n).exists()]
    if absent:
        raise FileNotFoundError(
            "Missing gridded CMIP6 array(s) in "
            f"{CMIP6}: {', '.join(absent)}.\n"
            "Download them from the figshare deposit (DOI 10.6084/m9.figshare.34419060, "
            "folder 07_CMIP6_gridded_fields/) and copy them into that directory. "
            "RUN.md lists the six files."
        )
    return [CMIP6 / n for n in names]
