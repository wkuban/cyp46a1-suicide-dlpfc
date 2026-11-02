# CYP46A1 in the prefrontal cortex of suicide decedents

Analysis code and results for:

> Kuban W, Pukło R, Giner S, García-Gutiérrez MS.
> *CYP46A1 expression in the dorsolateral prefrontal cortex of suicide decedents.*
> Manuscript in preparation (Journal of Psychiatric Research).

Everything reported in the paper can be regenerated from this repository: the
RT-qPCR statistics, the three public-dataset re-analyses, the normalisation
control, and every figure.

---

## What this study did

RT-qPCR of cytochrome P450 transcripts and two electron-transfer partners in the
dorsolateral prefrontal cortex (Brodmann area 9) of 12 male suicide decedents and
13 male natural-death controls, none with a documented psychiatric diagnosis,
psychotropic treatment or positive toxicology. Five targets were quantifiable
(CYP46A1, CYP2D6, CYP19A1, POR, FDX1); CYP11A1 and CYP17A1 amplified at or beyond
the limit of quantification and were not analysed.

CYP46A1 and CYP2D6, the two targets at the margin of nominal significance, were
then examined in three public datasets from the same brain region:

| Dataset | Design | What it adds |
|---|---|---|
| **GSE144136** | snRNA-seq, 17 male suicide + 17 male controls | cell-type resolution |
| **GSE101521** | bulk RNA-seq, 29 control / 9 MDD / 21 MDD-suicide | separates suicide from depression |
| **GSE213982** | snRNA-seq, 71 donors: 38 women plus the GSE144136 men, one pipeline | both sexes, largest cell-resolved comparison |

Three points of method that the code makes explicit:

1. **Every ΔCt is formed within one plate sample.** Script 10 rebuilds the RT-qPCR
   analysis directly from the QuantStudio exports, pairing each target Ct with
   the B2M Ct of the same sample. An earlier per-donor table had been assembled
   column by column and could pair one donor's target with another donor's B2M;
   that table and the statistics derived from it have been removed from this
   repository and are superseded by `results/qpcr_rebuilt_*`.
2. **All RT-qPCR tests run on the ΔCt (log₂) scale**, not on relative expression
   (2<sup>−ΔΔCt</sup>), which is right-skewed. No data points are excluded from
   the primary analysis.
3. **The GSE213982 matrix is not independent of GSE144136** — it contains the same
   male donors reprocessed. `data/metadata/gse213982_male_donor_matching.csv`
   documents the barcode-overlap matching that established this.

---

## Repository layout

```
code/
  02_gse144136_snrnaseq.py       GSE144136 re-analysis       -> Figure 6, Tables S5-S6
  03_gse101521_bulk.py           GSE101521 re-analysis       -> Figure 7, Tables S7-S8
  04_gse213982_pooled.py         GSE213982 re-analysis       -> Figure 8A-B, Tables S10-S11
  05_sokolov_lookup.py           meta-analysis table lookup
  06_normalisation_control.py    transcript-ratio control    -> Figure 8C, Table S12
  07_equivalence_and_power.py    TOST and power in GSE101521 -> Table S9
  08_cell_proportions.py         cell-class proportions      -> Figure 8D, Table S12
  09_qpcr_reference_gene_check.py  B2M and candidate reference genes from the exports -> Table S13
  10_qpcr_rebuild_from_exports.py  RT-qPCR statistics from the exports -> Table 2, Tables S2-S3, S13
  11_qpcr_sensitivity.py         Mann-Whitney and Grubbs sensitivity analyses -> Table S4
  figures/
    figures_2_to_5.py            Figures 2, 3, 4, 5
    figure_6.py  figure_7.py  figure_8.py

data/
  metadata/                      sample and donor tables (included)
  README.md                      how to obtain the public datasets (not included)

results/                         output of every script (included)
figures/                         generated figures (created on run)
```

## Running it

```bash
pip install -r requirements.txt

python code/11_qpcr_sensitivity.py       # runs from the committed per-donor table
python code/06_normalisation_control.py  # runs from the committed results
python code/07_equivalence_and_power.py

python code/figures/figures_2_to_5.py    # all figure scripts run from results/
python code/figures/figure_6.py
python code/figures/figure_7.py
python code/figures/figure_8.py
```

Scripts 02 to 05 and 08 need the public datasets first — see
**[data/README.md](data/README.md)**. Scripts 09 and 10 need the QuantStudio
`.xls` exports (`python code/10_qpcr_rebuild_from_exports.py data/qpcr_exports`),
which are not committed; their output is. All scripts resolve their own paths,
so they can be run from anywhere.

`04_gse213982_pooled.py` reads a 361-million-entry matrix and caches its
accumulators to `results/gse213982_accumulators.npz`, so the expensive pass
happens only once. Delete that file to force a full re-read.

## Data included here

`results/qpcr_rebuilt_per_donor.csv` contains, for every donor and gene, the
mean target Ct, the number and range of the replicates, the B2M Ct of the same
plate sample and the resulting ΔCt. Supplementary Table S2 (individual
2<sup>−ΔΔCt</sup> values) and all RT-qPCR statistics follow from this file, so
the group comparison is reproducible from the repository alone.

The instrument exports themselves are available from the corresponding author on
reasonable request.

## Software

Two environments were used, and the paper reports both:

| | Python | SciPy | NumPy | pandas |
|---|---|---|---|---|
| RT-qPCR statistics | 3.13.16 | 1.18.1 | 2.5.3 | — |
| Public-dataset analyses | 3.10.12 | 1.15.3 | 2.2.6 | 2.3.3 |

Results are not sensitive to these versions; `requirements.txt` pins only lower
bounds.

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
