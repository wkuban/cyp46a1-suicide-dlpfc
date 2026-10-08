#!/usr/bin/env python3
"""
Normalisation control for the pooled cell-resolved cohort (GSE213982 + the
GSE144136 male donors it contains).

An apparent difference in a single gene can arise from a shift in overall
transcript composition rather than from that gene. CYP46A1 is therefore
expressed per donor as a ratio to three neuronal reference transcripts measured
in the same nuclei (SNAP25, SLC17A7, RBFOX3), and the ratios are tested on the
log2 scale. Cell-class proportions are compared between groups for the same
reason.

These are single-gene control analyses and are reported without further
multiple-testing correction.

Input : results/gse213982_pseudobulk.csv
Output: results/gse213982_normalisation_control.csv
        results/gse213982_cell_proportions.csv

Reproduces Figure 8C-D and Supplementary Table S12.
"""
import os

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")

GENE = "CYP46A1"
REFERENCES = ["SNAP25", "SLC17A7", "RBFOX3"]
CLASSES = ["ExN", "InN", "ALL"]
PROPORTION_CLASSES = ["ExN", "InN", "Oli", "Ast", "OPC", "End", "Mic", "Mix"]


def welch_with_ci(a, b):
    """Welch t-test on log2 values; returns difference, 95% CI and p."""
    na, nb = len(a), len(b)
    d = b.mean() - a.mean()
    se = np.sqrt(a.var(ddof=1) / na + b.var(ddof=1) / nb)
    df = se**4 / ((a.var(ddof=1) / na) ** 2 / (na - 1)
                  + (b.var(ddof=1) / nb) ** 2 / (nb - 1))
    tc = stats.t.ppf(0.975, df)
    p = stats.ttest_ind(a, b, equal_var=False).pvalue
    return d, d - tc * se, d + tc * se, p


def main():
    P = pd.read_csv(os.path.join(RES, "gse213982_pseudobulk.csv"))

    # ---- ratios to neuronal reference transcripts ----
    rows = []
    for cls in CLASSES:
        W = (P[P.cell_class == cls]
             .pivot_table(index=["donor", "group"], columns="gene", values="log2CPM")
             .reset_index())
        for ref in [None] + REFERENCES:
            y = W[GENE] - (W[ref] if ref else 0)
            a = y[W.group == "Control"].values
            b = y[W.group == "Case"].values
            d, lo, hi, p = welch_with_ci(a, b)
            rows.append({
                "analysis": f"{GENE} / {ref}" if ref else f"{GENE}, no reference (CPM)",
                "cell_class": "All nuclei" if cls == "ALL" else cls,
                "n_control": len(a), "n_case": len(b),
                "ratio": 2**d, "ci_low": 2**lo, "ci_high": 2**hi, "p": p,
            })
    ratios = pd.DataFrame(rows)
    ratios.to_csv(os.path.join(RES, "gse213982_normalisation_control.csv"), index=False)

    # ---- cell-class proportions ----
    n = P[P.gene == GENE][["donor", "cell_class", "group", "n_nuclei"]]
    total = n[n.cell_class == "ALL"].set_index("donor").n_nuclei
    sub = n[n.cell_class != "ALL"].copy()
    sub["pct"] = sub.n_nuclei / sub.donor.map(total) * 100
    prop = []
    for cls in PROPORTION_CLASSES:
        a = sub.pct[(sub.cell_class == cls) & (sub.group == "Control")].values
        b = sub.pct[(sub.cell_class == cls) & (sub.group == "Case")].values
        prop.append({"cell_class": cls,
                     "control_pct": a.mean(), "case_pct": b.mean(),
                     "p": stats.ttest_ind(a, b, equal_var=False).pvalue})
    proportions = pd.DataFrame(prop)
    proportions.to_csv(os.path.join(RES, "gse213982_cell_proportions.csv"), index=False)

    pd.set_option("display.width", 200)
    print("Ratio to neuronal reference transcripts\n")
    print(ratios.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    print("\n\nCell-class proportions (% of nuclei per donor)\n")
    print(proportions.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    print(f"\nWritten to {RES}")


if __name__ == "__main__":
    main()
