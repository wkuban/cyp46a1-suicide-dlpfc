# Data

## Included in this repository

| Path | Contents |
|---|---|
| `qpcr/qpcr_relative_expression.csv` | Individual 2<sup>−ΔΔCt</sup> values, all seven targets, both groups. Same values as Supplementary Table S2. |
| `qpcr/prism_raw_values.json` | The same data as extracted from the original GraphPad Prism tables; read directly by the figure scripts. |
| `metadata/gse101521_samples.csv` | Diagnosis, age, sex, RIN, PMI and brain pH per sample, from the GSE101521 series matrix (platform GPL16791). |
| `metadata/gse213982_donors.csv` | Donor, group and sex for all 71 donors in the pooled cohort. |
| `metadata/gse213982_male_donor_matching.csv` | Barcode-overlap matching of the anonymous male labels in GSE213982 to their GSE144136 donor identity, with the fraction of shared barcodes. |

Everything in `results/` was produced from these files plus the public datasets below.

## Not included — download these first

The three expression datasets total roughly 1.5 GB and are not committed. Download
them from GEO into the folders shown. Scripts 02 to 05 expect exactly these paths.

### GSE144136 — snRNA-seq, male cohort (≈415 MB)

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144136

Into `data/gse144136/`, without unpacking:

```
GSE144136_GeneBarcodeMatrix_Annotated.mtx.gz
GSE144136_CellNames.csv.gz
GSE144136_GeneNames.csv.gz
```

Do not download `*premrna.tar.gz` (10 GB); it is not used.

### GSE101521 — bulk RNA-seq, three diagnostic groups (≈13 MB)

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE101521

Into `data/gse101521/`:

```
GSE101521_totalRNA_counts.csv.gz
```

Sample annotation is already provided in `metadata/gse101521_samples.csv`. It was
derived from `GSE101521-GPL16791_series_matrix.txt.gz`, available from the
`matrix/` subdirectory of the same accession if you wish to regenerate it. Note
that the counts matrix labels its columns `V1`–`V59`; these correspond to the
sample titles in the series matrix in the same order.

### GSE213982 — snRNA-seq, pooled female and male cohort (≈1.1 GB)

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213982

Into `data/gse213982/`, without unpacking:

```
GSE213982_combined_counts_matrix.mtx.gz
GSE213982_combined_counts_matrix_cells_columns.csv.gz
GSE213982_combined_counts_matrix_genes_rows.csv.gz
```

This matrix contains both the female donors and the male donors of GSE144136,
reprocessed through the same pipeline. The two analyses are therefore **not
independent**; see the manuscript and `metadata/gse213982_male_donor_matching.csv`.

### Sokolov et al. meta-analysis supplements (≈105 MB)

https://www.nature.com/articles/s41398-026-03978-8

Download the Supplementary Data files, unpack any archives, and place the
resulting tables in `data/sokolov/`. Script 05 scans every table it finds.

## A note on the MatrixMarket files

Both `.mtx.gz` matrices are 1-based. Scripts 02 and 04 decrement both the gene
and the cell index on reading. Getting this wrong shifts every gene by one row
and produces plausible-looking but entirely wrong results — the marker-gene
quality control in both scripts exists to catch exactly that.
