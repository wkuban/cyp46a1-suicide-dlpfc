#!/usr/bin/env python3
"""
Equivalence testing (TOST) and power for the bulk cohort (GSE101521).

The RT-qPCR study produced large standardised effects in a small cohort.
This script asks the complementary question of whether the much larger bulk
cohort is consistent with effects of that size. Two one-sided tests (TOST)
are run against equivalence bounds expressed on the standardised scale
(Lakens, 2017), and the power of the bulk comparison to detect effects of
those sizes is computed from the noncentral t distribution.

Bounds are the Hedges' g values observed by RT-qPCR for each gene
separately, plus a conventional medium bound of |g| = 0.5.

Input : results/gse101521_expression.csv   (committed)
Output: results/gse101521_equivalence_power.csv

Reproduces the TOST and power values reported in Results 3.7 and Table S9.
"""
import os

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")

# Hedges' g observed by RT-qPCR, per gene (see results/qpcr_statistics.json)
QPCR_G = {"CYP46A1": 1.08, "CYP2D6": 0.83}
EXTRA_BOUND = 0.5
ALPHA = 0.05


def welch(a, b):
    na, nb = len(a), len(b)
    va, vb = a.var(ddof=1), b.var(ddof=1)
    se = np.sqrt(va / na + vb / nb)
    diff = b.mean() - a.mean()
    df = (va / na + vb / nb) ** 2 / ((va / na) ** 2 / (na - 1) + (vb / nb) ** 2 / (nb - 1))
    sp = np.sqrt(((na - 1) * va + (nb - 1) * vb) / (na + nb - 2))
    return diff, se, df, sp, na, nb


def tost(diff, se, df, bound_raw):
    """Two one-sided tests. Returns (p_lower, p_upper, p_tost)."""
    p_lower = stats.t.sf((diff + bound_raw) / se, df)   # H0: diff <= -bound
    p_upper = stats.t.cdf((diff - bound_raw) / se, df)  # H0: diff >= +bound
    return p_lower, p_upper, max(p_lower, p_upper)


def power_two_sample(g, na, nb, alpha=ALPHA):
    """Power of a two-sided two-sample t test for standardised effect g."""
    df = na + nb - 2
    ncp = abs(g) * np.sqrt(na * nb / (na + nb))
    crit = stats.t.ppf(1 - alpha / 2, df)
    return (stats.nct.sf(crit, df, ncp) + stats.nct.cdf(-crit, df, ncp))


def main():
    E = pd.read_csv(os.path.join(RES, "gse101521_expression.csv"))
    rows = []
    for gene, g_qpcr in QPCR_G.items():
        d = E[E.gene == gene]
        a = d.log2CPM[d.diagnosis == "CON"].values
        b = d.log2CPM[d.diagnosis == "MDD-S"].values
        diff, se, df, sp, na, nb = welch(a, b)
        for label, g_bound in (("RT-qPCR estimate", g_qpcr), ("medium", EXTRA_BOUND)):
            bound_raw = g_bound * sp          # standardised bound -> log2 units
            pl, pu, p = tost(diff, se, df, bound_raw)
            rows.append({
                "gene": gene, "bound_type": label, "bound_g": g_bound,
                "bound_log2": bound_raw, "n_control": na, "n_case": nb,
                "log2_difference": diff, "fold_change": 2 ** diff,
                "se": se, "df": df,
                "p_lower": pl, "p_upper": pu, "p_TOST": p,
                "equivalent_at_0.05": p < ALPHA,
                "power_for_this_g": power_two_sample(g_bound, na, nb),
            })
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(RES, "gse101521_equivalence_power.csv"), index=False)

    pd.set_option("display.width", 200)
    show = ["gene", "bound_type", "bound_g", "n_control", "n_case",
            "fold_change", "p_TOST", "equivalent_at_0.05", "power_for_this_g"]
    print(out[show].to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    print(f"\nWritten to {os.path.join(RES, 'gse101521_equivalence_power.csv')}")


if __name__ == "__main__":
    main()
