#!/usr/bin/env python3
"""Rebuild the RT-qPCR analysis directly from the QuantStudio exports.

Why this exists: the per-donor tables used for the original analysis were
assembled one gene at a time from lists of available values, so for targets
with missing measurements a donor's target Ct could be paired with another
donor's B2M Ct. This script never leaves the donor: every Delta Ct is formed
from the target and B2M readings of the same plate sample, and a target
measurement with no B2M for that sample is dropped rather than filled.

Sample labels on the plates (C* control, S* suicide, NC no-template control)
carry group membership, which is all the group comparison needs; linking a
label to a case number in Table 1 is a separate record-keeping question and
does not enter the statistics.

Usage:  python 10_qpcr_rebuild_from_exports.py <dir with the .xls exports>

Output: results/qpcr_rebuilt_per_donor.csv   (per sample and gene)
        results/qpcr_rebuilt_summary.csv     (per gene, group comparison)
        results/qpcr_rebuilt.json
"""
import csv, glob, hashlib, json, os, re, sys, collections
import statistics as st

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import importlib.util
_spec = importlib.util.spec_from_file_location(
    "refcheck", os.path.join(HERE, "09_qpcr_reference_gene_check.py"))
refcheck = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(refcheck)

REFERENCE = "B2M"
# The seven targets prespecified in the manuscript. The group comparison and the
# Benjamini-Hochberg correction run over these only: CYC1, EIF4A2 and AKR1B1 were
# reference-gene candidates and further assays, not study targets, and belong in
# the assay-quality table instead.
PRESPECIFIED = ["CYP46A1", "CYP2D6", "CYP11A1", "CYP17A1", "CYP19A1", "POR", "FDX1"]
ARTEFACT_BELOW = 15.0          # Ct this low is an instrument artefact, not signal
# A target enters the group comparison only if it is actually quantifiable:
# Whether an assay is quantifiable is a property of the assay, not of how many
# plate exports happen to be available: a target is analysed if its signal is
# within the quantifiable range and its wells return a Ct, and the donors it
# covers are then reported as they are.
MAX_MEDIAN_CT = 35.0           # median Ct across donors
MIN_CALL_RATE = 0.80           # fraction of attempted wells returning a Ct
ASSAY = {  # label on the plate -> gene name
    "HKG1 B2M": "B2M", "HKG 3 CYC1": "CYC1", "HKG2 EIF4A2": "EIF4A2",
    "GOI 2 CYP46A1": "CYP46A1", "GOI CYP2D6": "CYP2D6", "GOI 2 CYP11A1": "CYP11A1",
    "GOI 2 CYP17A1": "CYP17A1", "GOI 3 CYP19A1": "CYP19A1", "GOI 9 POR": "POR",
    "GOI 10 FDX": "FDX1", "GOI 11 AKR1B1": "AKR1B1",
    "GOI 7 CYP11B2": "CYP11B2", "GOI 8 CYP21A2": "CYP21A2",
}


def read_exports(folder):
    """Every well from every export, duplicates of the same file read once."""
    seen, rows = set(), []
    for path in sorted(glob.glob(os.path.join(folder, "*.xls"))):
        digest = hashlib.md5(open(path, "rb").read()).hexdigest()
        if digest in seen:
            continue
        seen.add(digest)
        rows += refcheck.read_plate(path)
    return rows


def per_sample_ct(rows):
    """({gene: {sample: (mean Ct, n replicates, range)}}, {gene: call rate})"""
    wells = collections.defaultdict(lambda: collections.defaultdict(list))
    attempted = collections.Counter()
    called = collections.Counter()
    for r in rows:
        gene = ASSAY.get(r["target"])
        if gene is None or r["sample"] in ("", "NC"):
            continue
        attempted[gene] += 1
        if not re.match(r"^[\d.]+$", r["ct"] or ""):
            continue
        value = float(r["ct"])
        if value < ARTEFACT_BELOW:
            continue
        called[gene] += 1
        wells[gene][r["sample"]].append(value)
    out = {}
    for gene, samples in wells.items():
        out[gene] = {s: (st.mean(v), len(v), (max(v) - min(v)) if len(v) > 1 else 0.0)
                     for s, v in samples.items()}
    rate = {g: called[g] / attempted[g] for g in attempted if attempted[g]}
    return out, rate


def hedges_g(a, b):
    na, nb = len(a), len(b)
    sp = np.sqrt(((na - 1) * np.var(a, ddof=1) + (nb - 1) * np.var(b, ddof=1))
                 / (na + nb - 2))
    d = (np.mean(a) - np.mean(b)) / sp
    return d * (1 - 3 / (4 * (na + nb) - 9))


def benjamini_hochberg(pvals):
    p = np.asarray(pvals, float)
    m = len(p)
    order = np.argsort(p)
    q, running = np.empty(m), 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        running = min(running, p[i] * m / rank)
        q[i] = running
    return q


def compare(dct_control, dct_case):
    """Welch on the Delta Ct scale, with the fold change it implies."""
    a, b = list(dct_control.values()), list(dct_case.values())
    t, p = stats.ttest_ind(a, b, equal_var=False)
    try:
        p_mw = stats.mannwhitneyu(a, b).pvalue
    except ValueError:
        p_mw = float("nan")
    df = stats.ttest_ind(a, b, equal_var=False).df
    diff = np.mean(b) - np.mean(a)                 # Delta Delta Ct, case versus control
    se = np.sqrt(np.var(a, ddof=1) / len(a) + np.var(b, ddof=1) / len(b))
    crit = stats.t.ppf(0.975, df)
    return {
        "n_control": len(a), "n_case": len(b),
        "ddct": diff, "fold_change": 2 ** -diff,
        "fc_lo": 2 ** -(diff + crit * se), "fc_hi": 2 ** -(diff - crit * se),
        "t": abs(t), "df": df, "p": p, "p_mw": p_mw,
        "hedges_g": hedges_g(a, b), "eta2": t ** 2 / (t ** 2 + df),
        "sd_control": st.stdev(a), "sd_case": st.stdev(b),
    }


def main(folder):
    ct, call_rate = per_sample_ct(read_exports(folder))
    ref = ct.get(REFERENCE, {})
    if not ref:
        raise SystemExit("no B2M readings found in the exports")

    targets = [g for g in ct if g != REFERENCE]
    # quality is reported for every assay; only prespecified targets are tested
    per_donor, quality, stats_by_gene = [], {}, {}

    for gene in sorted(targets):
        dct = {"control": {}, "suicide": {}}
        for sample, (mean_ct, n_rep, spread) in sorted(ct[gene].items()):
            if sample not in ref:                   # no reference for this donor
                continue
            group = "control" if sample.upper().startswith("C") else "suicide"
            value = mean_ct - ref[sample][0]
            dct[group][sample] = value
            per_donor.append({
                "gene": gene, "sample": sample, "group": group,
                "ct_target": round(mean_ct, 3), "replicates": n_rep,
                "replicate_range": round(spread, 3),
                "ct_b2m": round(ref[sample][0], 3), "dct": round(value, 3),
            })
        med = st.median([v[0] for v in ct[gene].values()])
        donors = {g: len(d) for g, d in dct.items()}
        rate = call_rate.get(gene, 0.0)
        ok = med <= MAX_MEDIAN_CT and rate >= MIN_CALL_RATE
        quality[gene] = {"median_ct": round(med, 2), "call_rate": round(rate, 3),
                         "donors": donors, "prespecified": gene in PRESPECIFIED,
                         "analysed": bool(ok and gene in PRESPECIFIED)}
        if (ok and gene in PRESPECIFIED
                and donors["control"] >= 3 and donors["suicide"] >= 3):
            base = st.mean(dct["control"].values())
            res = compare(dct["control"], dct["suicide"])
            res["calibrator_dct"] = base
            stats_by_gene[gene] = res

    analysed = sorted(stats_by_gene)
    qs = benjamini_hochberg([stats_by_gene[g]["p"] for g in analysed])
    for g, q in zip(analysed, qs):
        stats_by_gene[g]["q_fdr"] = float(q)

    outdir = os.path.join(ROOT, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "qpcr_rebuilt_per_donor.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(per_donor[0].keys()))
        w.writeheader(); w.writerows(per_donor)
    with open(os.path.join(outdir, "qpcr_rebuilt_summary.csv"), "w", newline="") as f:
        cols = ["gene", "n_control", "n_case", "ddct", "fold_change", "fc_lo", "fc_hi",
                "t", "df", "p", "p_mw", "q_fdr", "hedges_g", "eta2"]
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for g in analysed:
            w.writerow({"gene": g, **{k: round(v, 6) if isinstance(v, float) else v
                                      for k, v in stats_by_gene[g].items()}})
    json.dump({"quality": quality, "statistics": stats_by_gene},
              open(os.path.join(outdir, "qpcr_rebuilt.json"), "w"), indent=1)

    print(f"{'gene':9} {'medCt':>6} {'call rate':>10} {'donors C/S':>11} {'analysed':>9}")
    for g in sorted(quality):
        q = quality[g]
        print(f"{g:9} {q['median_ct']:6.2f} {q['call_rate']:10.2f} "
              f"{q['donors']['control']:5d}/{q['donors']['suicide']:<5d} "
              f"{'yes' if q['analysed'] else ('no' if q['prespecified'] else 'not a target'):>12}")
    print(f"\nGroup comparison on Delta Ct, {len(analysed)} analysed targets")
    print(f"{'gene':9} {'n C/S':>7} {'ddCt':>7} {'fold change (95% CI)':>24} "
          f"{'t (df)':>13} {'p':>8} {'q':>7} {'g':>6}")
    for g in analysed:
        r = stats_by_gene[g]
        print(f"{g:9} {r['n_control']:3d}/{r['n_case']:<3d} {r['ddct']:7.2f} "
              f"{r['fold_change']:8.2f} ({r['fc_lo']:.2f}-{r['fc_hi']:.2f}){'':>4} "
              f"{r['t']:5.2f} ({r['df']:5.1f}) {r['p']:8.4f} {r['q_fdr']:7.4f} "
              f"{r['hedges_g']:6.2f}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "data/qpcr_exports")
