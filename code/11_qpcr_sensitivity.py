#!/usr/bin/env python3
"""Sensitivity analyses for the RT-qPCR group comparison (Table S4).

Reads the per-donor Delta Ct table written by 10_qpcr_rebuild_from_exports.py
and, for each analysed target, reports

  * the Mann-Whitney U test on Delta Ct (U given for the control group),
  * a two-sided Grubbs test (alpha = 0.05, single pass, within each group) and
    the Welch t-test after removing any flagged value.

Input : results/qpcr_rebuilt_per_donor.csv
Output: results/qpcr_sensitivity.csv
"""
import csv
import os

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RES = os.path.join(ROOT, "results")

TARGETS = ["CYP46A1", "CYP2D6", "CYP19A1", "POR", "FDX1"]
ALPHA = 0.05


def grubbs_index(x, alpha=ALPHA):
    """Index of the single value flagged by a two-sided Grubbs test, or None."""
    x = np.asarray(x, float)
    n = len(x)
    if n < 3:
        return None
    dev = np.abs(x - x.mean())
    g = dev.max() / x.std(ddof=1)
    t = stats.t.ppf(1 - alpha / (2 * n), n - 2)
    g_crit = (n - 1) / np.sqrt(n) * np.sqrt(t ** 2 / (n - 2 + t ** 2))
    return int(dev.argmax()) if g > g_crit else None


def main():
    rows = [r for r in csv.DictReader(open(os.path.join(RES, "qpcr_rebuilt_per_donor.csv")))
            if r["dct"]]
    out = []
    for gene in TARGETS:
        ctrl = [(r["sample"], float(r["dct"])) for r in rows if r["gene"] == gene and r["group"] == "control"]
        sui = [(r["sample"], float(r["dct"])) for r in rows if r["gene"] == gene and r["group"] == "suicide"]
        c = [v for _, v in ctrl]
        s = [v for _, v in sui]
        p_welch = stats.ttest_ind(c, s, equal_var=False).pvalue
        mw = stats.mannwhitneyu(c, s, alternative="two-sided")
        ic, is_ = grubbs_index(c), grubbs_index(s)
        flagged = []
        if ic is not None:
            flagged.append(f"{ctrl[ic][0]} (control)")
        if is_ is not None:
            flagged.append(f"{sui[is_][0]} (suicide)")
        c2 = [v for i, v in enumerate(c) if i != ic]
        s2 = [v for i, v in enumerate(s) if i != is_]
        p_after = stats.ttest_ind(c2, s2, equal_var=False).pvalue
        fc_after = 2 ** -(np.mean(s2) - np.mean(c2))
        out.append({
            "gene": gene,
            "n_control": len(c), "n_suicide": len(s),
            "p_welch_all": round(p_welch, 4),
            "U_control": int(mw.statistic), "p_mannwhitney": round(mw.pvalue, 4),
            "grubbs_flagged": "; ".join(flagged) if flagged else "none",
            "n_control_after_grubbs": len(c2), "n_suicide_after_grubbs": len(s2),
            "p_welch_after_grubbs": round(p_after, 4),
            "fold_change_after_grubbs": round(fc_after, 3),
        })
    path = os.path.join(RES, "qpcr_sensitivity.csv")
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    for r in out:
        print(f"{r['gene']:8s} U={r['U_control']:3d} p_MW={r['p_mannwhitney']:.3f}  "
              f"Grubbs: {r['grubbs_flagged']:22s} n={r['n_control_after_grubbs']}/{r['n_suicide_after_grubbs']} "
              f"p={r['p_welch_after_grubbs']:.3f}")
    print("wrote", path)


if __name__ == "__main__":
    main()
