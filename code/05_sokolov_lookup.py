#!/usr/bin/env python3
"""
Lookup of CYP46A1 and CYP2D6 in the supplementary tables of a 16-cohort
meta-analysis of post-mortem suicide brain transcriptomes
(Sokolov et al., Transl Psychiatry 2026;16(1):222, doi:10.1038/s41398-026-03978-8).

Those tables list only genes reaching significance, not genome-wide summary
statistics, so absence of a gene means it did not reach the reporting threshold
rather than that it was unchanged. Five genes the study itself reports as
significant are used as positive controls to confirm that the relevant tables
were read correctly.

Download the Supplementary Data files from the article page, unpack any archives,
and place them in data/sokolov/.

Input : data/sokolov/
Output: results/sokolov_hits.csv
        results/sokolov_table_inventory.csv
"""
import argparse
import glob
import os

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTDIR = os.path.join(ROOT, "results")

# Genes this meta-analysis itself reports as significant; if none of these are
# found either, the wrong files have been read.
POSITIVE_CONTROLS = ["P2RY12", "CX3CR1", "GPR34", "SOX9", "PMP2"]

SYMBOL_PATTERN = r"[A-Z][A-Z0-9\-\.]{1,19}"


def iter_tables(folder):
    """Yield (path, sheet_name, DataFrame) for every readable table in folder."""
    exts = (".xlsx", ".xls", ".xlsm", ".csv", ".tsv", ".txt")
    files = [f for f in glob.glob(os.path.join(folder, "**", "*"), recursive=True)
             if os.path.isfile(f) and f.lower().endswith(exts)
             and not os.path.basename(f).startswith("~$")]
    if not files:
        print("No .xlsx/.csv/.tsv/.txt file found in that folder. Contents:")
        for f in sorted(glob.glob(os.path.join(folder, "**", "*"), recursive=True))[:40]:
            print("   ", os.path.relpath(f, folder))
    for f in sorted(files):
        ext = os.path.splitext(f)[1].lower()
        try:
            if ext in (".xlsx", ".xls", ".xlsm"):
                xl = pd.ExcelFile(f)
                for sh in xl.sheet_names:
                    try:
                        yield f, sh, xl.parse(sh, dtype=str)
                    except Exception as e:
                        print(f"  (could not read {os.path.basename(f)} / {sh}: {e})")
            else:
                sep = "\t" if ext in (".tsv", ".txt") else ","
                yield f, "-", pd.read_csv(f, sep=sep, dtype=str, engine="python",
                                          on_bad_lines="skip")
        except Exception as e:
            print(f"  (could not read {os.path.basename(f)}: {e})")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dir", default=os.path.join(ROOT, "data", "sokolov"),
                    help="folder containing the downloaded supplementary tables")
    ap.add_argument("--genes", nargs="+", default=["CYP46A1", "CYP2D6"])
    args = ap.parse_args()
    os.makedirs(OUTDIR, exist_ok=True)

    if not os.path.isdir(args.dir) or not os.listdir(args.dir):
        os.makedirs(args.dir, exist_ok=True)
        print("\n" + "=" * 72)
        print("The data folder is empty. It has just been created at:")
        print("   ", os.path.abspath(args.dir))
        print("\nDownload the Supplementary Data files from")
        print("    https://www.nature.com/articles/s41398-026-03978-8")
        print("unpack any archives, place the files in that folder, and re-run.")
        print("=" * 72)
        return

    targets = {g.upper() for g in args.genes}
    controls = {g.upper() for g in POSITIVE_CONTROLS}
    hits, inventory = [], []
    all_symbols, found_controls = set(), set()
    n_tables = 0

    print(f"Folder: {os.path.abspath(args.dir)}\n")
    for path, sheet, df in iter_tables(args.dir):
        n_tables += 1
        name = os.path.relpath(path, args.dir)
        if df is None or df.empty:
            inventory.append({"file": name, "sheet": sheet, "rows": 0, "columns": 0})
            print(f"  [empty]  {name} / {sheet}")
            continue
        df.columns = [str(c) for c in df.columns]
        inventory.append({"file": name, "sheet": sheet, "rows": len(df),
                          "columns": len(df.columns),
                          "headers": " | ".join(df.columns[:8])})
        print(f"  {len(df):>7} x {len(df.columns):<3}  {name} / {sheet}")

        # Search every column: gene symbols are not consistently placed
        for col in df.columns:
            vals = df[col].astype(str).str.strip().str.upper()
            symbols = vals[vals.str.fullmatch(SYMBOL_PATTERN) == True]  # noqa: E712
            all_symbols.update(symbols.unique()[:200000])
            found_controls.update(set(symbols.unique()) & controls)
            mask = vals.isin(targets)
            if mask.any():
                for _, r in df[mask].iterrows():
                    rec = {"file": name, "sheet": sheet, "column": col,
                           "gene": str(r[col]).strip().upper()}
                    for c in df.columns:
                        if c != col:
                            rec[c] = r[c]
                    hits.append(rec)

    pd.DataFrame(inventory).to_csv(
        os.path.join(OUTDIR, "sokolov_table_inventory.csv"),
        index=False, encoding="utf-8-sig")
    print(f"\nTables examined: {n_tables}")

    print("\n" + "-" * 72)
    print("POSITIVE CONTROLS (genes this study reports as significant):")
    for g in POSITIVE_CONTROLS:
        print(f"   {'FOUND  ' if g in found_controls else 'missing'}  {g}")
    print("-" * 72)

    if hits:
        out = pd.DataFrame(hits)
        out.to_csv(os.path.join(OUTDIR, "sokolov_hits.csv"),
                   index=False, encoding="utf-8-sig")
        print(f"\n{len(out)} rows containing the target genes. "
              f"Written to results/sokolov_hits.csv\n")
        with pd.option_context("display.width", 220, "display.max_columns", 40):
            print(out.head(30).to_string(index=False))
        return

    print("\nNot found:", ", ".join(sorted(targets)))
    if n_tables == 0:
        print("No table could be read at all - the files are probably missing "
              "or still inside an archive.")
    elif not found_controls:
        print("None of the positive-control genes was found either, so these files "
              "do not contain gene lists. Download the remaining Supplementary Data "
              "files and re-run.")
    else:
        print(f"Positive controls found ({', '.join(sorted(found_controls))}), so the "
              f"gene tables were read correctly, but neither target gene appears in them.")
        print(f"Unique gene symbols across all tables: approximately {len(all_symbols)}.")
        if len(all_symbols) > 10000:
            print("That looks like a genome-wide list, so absence would be informative.")
        else:
            print("That looks like a significant-genes-only list, so absence means the "
                  "genes did not reach the reporting threshold, not that they were "
                  "unchanged.")


if __name__ == "__main__":
    main()
