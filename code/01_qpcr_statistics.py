#!/usr/bin/env python3
"""
RT-qPCR statistics: cytochrome P450 expression in human DLPFC (BA9),
suicide decedents versus natural-death controls.

All tests are performed on the DeltaCt (log2) scale rather than on relative
expression (2^-DDCt), which is right-skewed. The input file contains the
relative expression values exactly as reported in Supplementary Table S2;
taking log2 of these values recovers the -DDCt scale used for inference.

No data points are excluded from the primary analysis. A Grubbs-based
sensitivity analysis is reported separately.

Input : data/qpcr/qpcr_relative_expression.csv
Output: results/qpcr_statistics.json
        results/qpcr_sensitivity.json

Reproduces Table 2 and Supplementary Tables S3 and S4 of the manuscript.
"""
import csv
import json
import os
from collections import defaultdict

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
IN = os.path.join(ROOT, "data", "qpcr", "qpcr_relative_expression.csv")
OUTDIR = os.path.join(ROOT, "results")

# Order used for Benjamini-Hochberg correction across the seven prespecified targets
GENES = ["CYP2D6", "CYP46A1", "CYP11A1", "CYP17A1", "CYP19A1", "POR", "FDX1"]


def load():
    d = defaultdict(lambda: {"control": [], "suicide": []})
    with open(IN, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            d[r["gene"]][r["group"]].append(float(r["relative_expression"]))
    return d


def benjamini_hochberg(p):
    p = np.asarray(p, float)
    m = len(p)
    order = np.argsort(p)
    q = np.empty(m)
    running = 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        running = min(running, p[i] * m / rank)
        q[i] = running
    return q


def grubbs(x, alpha=0.05):
    """Two-sided Grubbs test, single pass. Returns (index, G, G_critical)."""
    x = np.asarray(x, float)
    n = len(x)
    if n < 3:
        return None, np.nan, np.nan
    sd = x.std(ddof=1)
    if sd == 0:
        return None, np.nan, np.nan
    dev = np.abs(x - x.mean())
    i = int(np.argmax(dev))
    g = dev[i] / sd
    t = stats.t.ppf(1 - alpha / (2 * n), n - 2)
    g_crit = (n - 1) / np.sqrt(n) * np.sqrt(t**2 / (n - 2 + t**2))
    return i, float(g), float(g_crit)


def welch(c, s):
    """Welch t-test on the log2 scale; returns the full set of reported statistics."""
    nc, ns = len(c), len(s)
    vc, vs = c.var(ddof=1), s.var(ddof=1)
    se = np.sqrt(vc / nc + vs / ns)
    diff = s.mean() - c.mean()                      # log2 fold change, suicide vs control
    df = (vc / nc + vs / ns) ** 2 / ((vc / nc) ** 2 / (nc - 1) + (vs / ns) ** 2 / (ns - 1))
    t = diff / se
    p = 2 * stats.t.sf(abs(t), df)
    tcrit = stats.t.ppf(0.975, df)
    lo, hi = diff - tcrit * se, diff + tcrit * se
    eta2 = t * t / (t * t + df)
    # Hedges' g, with the small-sample correction J
    sp = np.sqrt(((nc - 1) * vc + (ns - 1) * vs) / (nc + ns - 2))
    g = (diff / sp) * (1 - 3 / (4 * (nc + ns) - 9))
    return dict(nc=nc, ns=ns, diff=diff, se=se, df=df, t=t, p=p,
                fc=2**diff, fc_lo=2**lo, fc_hi=2**hi, eta2=eta2, g=g, sp=sp)


def geometric_mean_ci(x):
    """Back-transformed mean and 95% CI of log2 values."""
    m = x.mean()
    h = stats.t.ppf(0.975, len(x) - 1) * x.std(ddof=1) / np.sqrt(len(x))
    return 2**m, 2 ** (m - h), 2 ** (m + h)


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    data = load()
    primary, sens = {}, {}

    for gene in GENES:
        # log2 of relative expression == -DDCt
        c = np.log2(np.array(data[gene]["control"], float))
        s = np.log2(np.array(data[gene]["suicide"], float))
        w = welch(c, s)

        primary[gene] = dict(
            nc=w["nc"], ns=w["ns"],
            # DDCt means; positive = lower expression in that group
            ddct_control_mean=-c.mean(), ddct_control_sd=c.std(ddof=1),
            ddct_suicide_mean=-s.mean(), ddct_suicide_sd=s.std(ddof=1),
            geomean_ci_control=geometric_mean_ci(c),
            geomean_ci_suicide=geometric_mean_ci(s),
            dct_difference=-w["diff"],
            t=abs(w["t"]), df=w["df"], p=w["p"],
            fold_change=w["fc"], fold_change_lo=w["fc_lo"], fold_change_hi=w["fc_hi"],
            eta2=w["eta2"], hedges_g=w["g"],
            shapiro_control=stats.shapiro(c).pvalue,
            shapiro_suicide=stats.shapiro(s).pvalue,
            levene=stats.levene(c, s).pvalue,
            mannwhitney_p=stats.mannwhitneyu(c, s).pvalue,
        )

        # Sensitivity: Welch test after removing any value flagged by Grubbs
        ic, gc, gcrit_c = grubbs(c)
        iss, gs, gcrit_s = grubbs(s)
        c2 = np.delete(c, ic) if (ic is not None and gc > gcrit_c) else c
        s2 = np.delete(s, iss) if (iss is not None and gs > gcrit_s) else s
        mw = stats.mannwhitneyu(c, s)
        sens[gene] = dict(
            grubbs_control=[ic, gc, gcrit_c], grubbs_suicide=[iss, gs, gcrit_s],
            n_after=[len(c2), len(s2)],
            p_after_grubbs=stats.ttest_ind(c2, s2, equal_var=False).pvalue,
            fold_change_after_grubbs=2 ** (s2.mean() - c2.mean()),
            mannwhitney_p=mw.pvalue, mannwhitney_U=mw.statistic,
        )

    q = benjamini_hochberg([primary[g]["p"] for g in GENES])
    for g, qq in zip(GENES, q):
        primary[g]["q_fdr"] = float(qq)

    json.dump(primary, open(os.path.join(OUTDIR, "qpcr_statistics.json"), "w"),
              indent=1, default=float)
    json.dump(sens, open(os.path.join(OUTDIR, "qpcr_sensitivity.json"), "w"),
              indent=1, default=float)

    print(f"{'gene':9s} {'n C/S':9s} {'FC':>6s} {'95% CI':>15s} "
          f"{'p':>8s} {'q':>7s} {'eta2':>6s} {'g':>7s} {'MW p':>8s}")
    print("-" * 80)
    for g in GENES:
        r = primary[g]
        ci = f"({r['fold_change_lo']:.2f}-{r['fold_change_hi']:.2f})"
        print(f"{g:9s} {r['nc']:2d}/{r['ns']:<6d} {r['fold_change']:6.3f} {ci:>15s} "
              f"{r['p']:8.4f} {r['q_fdr']:7.3f} {r['eta2']:6.3f} "
              f"{r['hedges_g']:7.3f} {r['mannwhitney_p']:8.4f}")
    print(f"\nWritten to {OUTDIR}")


if __name__ == "__main__":
    main()
