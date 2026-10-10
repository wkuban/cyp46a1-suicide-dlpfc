# CYP46A1 in the prefrontal cortex of suicide decedents

Analysis code and results for:

> Kuban W, Pukło R, Giner S, García-Gutiérrez MS.
> *CYP46A1 expression in the prefrontal cortex of suicide decedents and public
> transcriptomic cohorts.*
> Manuscript in preparation.

Everything reported in the paper can be regenerated from this repository: the
RT-qPCR statistics, the three public-dataset re-analyses, the normalisation
control, and every figure.

---

## What this study did

RT-qPCR of five cytochrome P450 enzymes and two electron-transfer partners in the
dorsolateral prefrontal cortex (Brodmann area 9) of 12 male suicide decedents and
13 male natural-death controls without psychiatric diagnosis, psychotropic
treatment or positive toxicology.

The two nominally positive findings were then tested against three public
datasets from the same brain region:

| Dataset | Design | What it adds |
|---|---|---|
| **GSE144136** | snRNA-seq, 17 male suicide + 17 male controls | cell-type resolution |
| **GSE101521** | bulk RNA-seq, 29 control / 9 MDD / 21 MDD-suicide | separates suicide from depression |
| **GSE213982** | snRNA-seq, 71 donors: 38 women plus the GSE144136 men, one pipeline | both sexes, largest cell-resolved comparison |

Two points of method that the code makes explicit:

1. **All RT-qPCR tests run on the ΔCt (log₂) scale**, not on relative expression
   (2<sup>−ΔΔCt</sup>), which is right-skewed. No data points are excluded from
   the primary analysis.
2. **The GSE213982 matrix is not independent of GSE144136** — it contains the same
   male donors reprocessed. `data/metadata/gse213982_male_donor_matching.csv`
   documents the barcode-overlap matching that established this.

---

## Repository layout

```
code/
  01_qpcr_statistics.py          RT-qPCR statistics        -> Table 2, Tables S3-S4
  02_gse144136_snrnaseq.py       GSE144136 re-analysis     -> Figure 6, Tables S5-S6
  03_gse101521_bulk.py           GSE101521 re-analysis     -> Figure 7, Tables S7-S9
  04_gse213982_pooled.py         GSE213982 re-analysis     -> Figure 8A-B, Tables S10-S11
  05_sokolov_lookup.py           meta-analysis table lookup
  06_normalisation_control.py    transcript-ratio control  -> Figure 8C, Table S12
  figures/
    figures_2_to_4.py            Figures 2, 3, 4, 5
    figure_6.py  figure_7.py  figure_8.py

data/
  qpcr/                          RT-qPCR values (included)
  metadata/                      sample and donor tables (included)
  README.md                      how to obtain the public datasets (not included)

results/                         output of every script (included)
figures/                         generated figures (created on run)
```

## Running it

```bash
pip install -r requirements.txt

python code/01_qpcr_statistics.py        # runs immediately, no download needed
python code/06_normalisation_control.py  # runs from the committed results

python code/figures/figures_2_to_4.py
python code/figures/figure_7.py
python code/figures/figure_8.py
```

Scripts 02 to 05 need the public datasets first — see **[data/README.md](data/README.md)**.
They resolve their own paths, so they can be run from anywhere.

`04_gse213982_pooled.py` reads a 361-million-entry matrix and caches its
accumulators to `results/gse213982_accumulators.npz`, so the expensive pass
happens only once. Delete that file to force a full re-read.

## Data included here

`data/qpcr/qpcr_relative_expression.csv` contains the individual
2<sup>−ΔΔCt</sup> values for all seven targets in both groups — the same values
published as Supplementary Table S2. Taking log₂ of these recovers the −ΔΔCt
scale used for inference, so the entire RT-qPCR analysis is reproducible from
this repository alone.

`data/qpcr/prism_raw_values.json` is the same data as extracted from the original
GraphPad Prism tables, kept because the figure scripts read it directly.

Raw Ct values are available from the corresponding author on reasonable request.

## Software

Two environments were used, and the paper reports both:

| | Python | SciPy | NumPy | pandas |
|---|---|---|---|---|
| RT-qPCR statistics | 3.13.16 | 1.18.1 | 2.5.3 | — |
| Public-dataset analyses | 3.10.12 | 1.15.3 | 2.2.6 | 2.3.3 |

The RT-qPCR results were additionally verified in GraphPad Prism 10. Results are
not sensitive to these versions; `requirements.txt` pins only lower bounds.

## Citing

Please cite the paper. If you use the code itself, see `CITATION.cff`.

## Licence

Code is released under the MIT Licence (`LICENSE`). The public datasets
re-analysed here remain under the terms set by their original depositors and
must be cited directly:

- Nagy C, et al. *Nat Neurosci* 2020;23(6):771-781. doi:10.1038/s41593-020-0621-y — GSE144136
- Pantazatos SP, et al. *Mol Psychiatry* 2017;22(5):760-773. doi:10.1038/mp.2016.130 — GSE101521
- Maitra M, et al. *Nat Commun* 2023;14(1):2912. doi:10.1038/s41467-023-38530-5 — GSE213982
- Sokolov AV, et al. *Transl Psychiatry* 2026;16(1):222. doi:10.1038/s41398-026-03978-8
