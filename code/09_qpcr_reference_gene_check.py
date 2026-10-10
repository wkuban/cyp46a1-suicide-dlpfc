"""Reference-gene assessment and B2M stability from QuantStudio 12K Flex exports.

Reads the instrument's own .xls exports (Results sheet: Well Position, Sample Name,
Target Name, CT) and reproduces Supplementary Table S13.

Usage:  python 09_qpcr_reference_gene_check.py <dir with *.xls>

Sample names follow the plate convention: C* = control donor, S* = suicide decedent,
NC = no-template control. Wells left unnamed in the plate setup are excluded.
Ct below 15 cycles is treated as an instrument artefact, not signal.
"""
import sys, os, re, glob, collections, statistics as st
import xlrd
from scipy import stats

ARTEFACT_BELOW = 15.0


def read_plate(path):
    """Return (rows, meta). Header row is located by name, not position, because
    the exports differ in how much metadata precedes it."""
    book = xlrd.open_workbook(path)
    sheet = book.sheet_by_name("Results")
    head = None
    for r in range(sheet.nrows):
        vals = [str(sheet.cell_value(r, c)).strip() for c in range(sheet.ncols)]
        if "Well Position" in vals:
            head, hdr = r, vals
            break
    if head is None:
        raise ValueError(f"no header row in {path}")
    idx = {name: i for i, name in enumerate(hdr) if name}
    rows = []
    for r in range(head + 1, sheet.nrows):
        v = [str(sheet.cell_value(r, c)).strip() for c in range(sheet.ncols)]
        well = v[idx["Well Position"]]
        if not re.match(r"^[A-H]\d{1,2}$", well or ""):
            continue
        rows.append({"plate": os.path.basename(path),
                     "well": well,
                     "sample": v[idx["Sample Name"]],
                     "target": v[idx["Target Name"]],
                     "ct": v[idx["CT"]]})
    return rows


def numeric(ct):
    return bool(re.match(r"^[\d.]+$", ct or ""))


def per_assay(rows):
    """Group the named-sample wells of one assay by donor, keeping the counts
    that describe data quality: dropouts, artefacts, replicate spread."""
    by_donor = collections.defaultdict(list)
    undetermined = artefacts = 0
    for r in rows:
        if r["sample"] in ("", "NC"):
            continue
        if not numeric(r["ct"]):
            undetermined += 1
            continue
        value = float(r["ct"])
        if value < ARTEFACT_BELOW:
            artefacts += 1
            continue
        by_donor[r["sample"]].append(value)
    spreads = [max(v) - min(v) for v in by_donor.values() if len(v) > 1]
    flat = [x for v in by_donor.values() for x in v]
    return {"donors": by_donor, "undetermined": undetermined, "artefacts": artefacts,
            "wells": len(flat), "median_ct": st.median(flat) if flat else None,
            "spread_median": st.median(spreads) if spreads else None,
            "spread_over_1": sum(1 for d in spreads if d > 1),
            "spread_n": len(spreads)}


def group_stability(by_donor):
    """Welch, Mann-Whitney and Levene on per-donor mean Ct, control vs suicide."""
    means = {k: st.mean(v) for k, v in by_donor.items()}
    ctl = [m for k, m in means.items() if k.upper().startswith("C")]
    case = [m for k, m in means.items() if k.upper().startswith("S")]
    if len(ctl) < 3 or len(case) < 3:
        return None
    out = {"n_control": len(ctl), "n_case": len(case),
           "mean_control": st.mean(ctl), "sd_control": st.stdev(ctl),
           "mean_case": st.mean(case), "sd_case": st.stdev(case),
           "p_welch": stats.ttest_ind(ctl, case, equal_var=False).pvalue,
           "p_mw": stats.mannwhitneyu(ctl, case).pvalue,
           "p_levene": stats.levene(ctl, case).pvalue}
    out["cv_control"] = 100 * out["sd_control"] / out["mean_control"]
    out["cv_case"] = 100 * out["sd_case"] / out["mean_case"]
    return out


def main(folder):
    rows = []
    for path in sorted(glob.glob(os.path.join(folder, "*.xls"))):
        rows += read_plate(path)
    if not rows:
        raise SystemExit(f"no .xls exports found in {folder}")

    assays = collections.defaultdict(list)
    for r in rows:
        if r["target"]:
            assays[r["target"]].append(r)

    print(f"{'assay':16} {'wells':>5} {'medCt':>6} {'undet':>6} {'artf':>5} "
          f"{'spread':>7} {'>1cyc':>7}  NTC")
    for name in sorted(assays):
        a = per_assay(assays[name])
        ntc = [r["ct"] for r in assays[name] if r["sample"] == "NC"]
        ntc = (", ".join(f"{float(x):.2f}" if numeric(x) else "undet" for x in ntc)
               or "none on plate")
        med = f"{a['median_ct']:.2f}" if a["median_ct"] else "-"
        spr = f"{a['spread_median']:.2f}" if a["spread_median"] else "-"
        print(f"{name:16} {a['wells']:5d} {med:>6} {a['undetermined']:6d} "
              f"{a['artefacts']:5d} {spr:>7} {a['spread_over_1']:3d}/{a['spread_n']:<3} {ntc}")

    print("\nGroup stability (per-donor mean Ct, control vs suicide):")
    for name in sorted(assays):
        s = group_stability(per_assay(assays[name])["donors"])
        if not s:
            continue
        print(f"  {name}: control {s['mean_control']:.2f} +/- {s['sd_control']:.2f} "
              f"(n={s['n_control']}, CV {s['cv_control']:.1f}%), "
              f"suicide {s['mean_case']:.2f} +/- {s['sd_case']:.2f} "
              f"(n={s['n_case']}, CV {s['cv_case']:.1f}%)")
        print(f"      Welch p = {s['p_welch']:.3f}, Mann-Whitney p = {s['p_mw']:.3f}, "
              f"Levene p = {s['p_levene']:.3f}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "data/qpcr_exports")
