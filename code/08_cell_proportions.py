#!/usr/bin/env python3
"""
Cell-class proportions in the pooled cohort (GSE213982), computed from the
complete nucleus annotation.

The pseudobulk table written by 04_gse213982_pooled.py keeps a donor x class
combination only when that donor contributes at least 20 nuclei of the class,
because fewer nuclei give an unstable expression estimate. That threshold is
appropriate for expression but not for composition: a donor with few nuclei of
a class has a low proportion of it, and dropping such donors truncates the
distribution from below and removes unequal numbers of donors per group.

Proportions are therefore computed here directly from the cell annotation,
before any expression filter, so that every donor contributes to every class.

Input : data/gse213982/GSE213982_combined_counts_matrix_cells_columns.csv.gz
        data/metadata/gse213982_donors.csv
Output: results/gse213982_nuclei_counts.csv      (per donor x class, unfiltered)
        results/gse213982_cell_proportions.csv   (group comparison)

Reproduces Figure 8D and the cell-proportion rows of Supplementary Table S12.
"""
import argparse
import gzip
import os
from collections import Counter

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
DEFAULT_CELLS = os.path.join(ROOT, "data", "gse213982",
                             "GSE213982_combined_counts_matrix_cells_columns.csv.gz")
DEFAULT_META = os.path.join(ROOT, "data", "metadata", "gse213982_donors.csv")


def parse(name):
    """Nucleus identifier -> (donor, broad cell class). Same rule as script 04."""
    p = name.split(".")
    d = p[0].replace("M24_2", "M24")
    return (d, p[2]) if len(p) >= 4 else (d, "?")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cells", default=DEFAULT_CELLS)
    ap.add_argument("--meta", default=DEFAULT_META)
    args = ap.parse_args()
    os.makedirs(RES, exist_ok=True)

    with gzip.open(args.cells, "rt") as fh:
        cells = [l.strip().strip('"') for l in fh][1:]
    meta = pd.read_csv(args.meta)
    meta["donor"] = meta.donor.str.replace("M24_2", "M24", regex=False)
    meta = meta.drop_duplicates("donor")
    grp = dict(zip(meta.donor, meta.group))
    sex = dict(zip(meta.donor, meta.sex))

    counts = Counter()
    total = Counter()
    for c in cells:
        d, b = parse(c)
        counts[(d, b)] += 1
        total[d] += 1
    donors = sorted(total)
    classes = sorted({b for _, b in counts})
    missing = [d for d in donors if d not in grp]
    if missing:
        raise SystemExit(f"donors missing from metadata: {missing}")

    rows = []
    for d in donors:
        for b in classes:
            n = counts.get((d, b), 0)           # zero is a real value, not a gap
            rows.append({"donor": d, "group": grp[d], "sex": sex[d],
                         "cell_class": b, "n_nuclei": n,
                         "n_nuclei_total": total[d],
                         "pct": n / total[d] * 100})
    N = pd.DataFrame(rows)
    N.to_csv(os.path.join(RES, "gse213982_nuclei_counts.csv"), index=False)

    out = []
    for b in classes:
        s = N[N.cell_class == b]
        a = s.pct[s.group == "Control"].values
        c = s.pct[s.group == "Case"].values
        t = stats.ttest_ind(a, c, equal_var=False)
        u = stats.mannwhitneyu(a, c)
        out.append({"cell_class": b, "n_control": len(a), "n_case": len(c),
                    "control_pct": a.mean(), "case_pct": c.mean(),
                    "control_sd": a.std(ddof=1), "case_sd": c.std(ddof=1),
                    "donors_with_zero": int((s.n_nuclei == 0).sum()),
                    "donors_under_20": int((s.n_nuclei < 20).sum()),
                    "p_welch": t.pvalue, "p_MW": u.pvalue})
    P = pd.DataFrame(out).sort_values("control_pct", ascending=False)
    P.to_csv(os.path.join(RES, "gse213982_cell_proportions.csv"), index=False)

    pd.set_option("display.width", 200)
    print(f"nuclei {len(cells)}, donors {len(donors)} "
          f"({sum(grp[d]=='Control' for d in donors)} control / "
          f"{sum(grp[d]=='Case' for d in donors)} case), classes {len(classes)}\n")
    print(P.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    print(f"\nWritten to {RES}")


if __name__ == "__main__":
    main()
