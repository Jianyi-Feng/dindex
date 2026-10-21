"""dindex: a bounded two-group balance index for microbial community composition.

    D = (P - M) / (P + M) = tanh(0.5 ln(P/M))

Quick start
-----------
>>> import pandas as pd, dindex
>>> counts = pd.DataFrame({"otu1": [10, 0], "otu2": [5, 20], "otu3": [1, 1]}, index=["s1", "s2"])
>>> taxonomy = {"otu1": "Imitervirales", "otu2": "Algavirales", "otu3": "Pimascovirales"}
>>> dindex.d_from_table(counts, taxonomy, plus="Imitervirales", minus="Algavirales")[["D", "share", "qc"]].round(3)
        D  share  qc
s1  0.333  0.938  ok
s2 -1.000  0.952  ok
"""
from .core import aggregate_groups, d_ci_bootstrap, d_ci_counts, d_from_log_ratio, d_from_table, d_index, log_ratio

__version__ = "0.1.0"
__all__ = ["d_index", "log_ratio", "d_from_log_ratio", "aggregate_groups", "d_from_table", "d_ci_counts", "d_ci_bootstrap", "__version__"]
