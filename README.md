# ESCC-STR-xQTL

Figure reproduction code for the multi-omics STR-xQTL study of esophageal squamous cell carcinoma.


## Scope

This repository contains custom R scripts for reproducing manuscript figures from processed analysis outputs. It is intended to support code review and figure reproduction for the manuscript.


## Repository contents

```text
.
|-- README.md
|-- ESCC-STR-xQTL_code.R
`-- outputs/                        # created when running the scripts
```

Current main script:

```text
ESCC-STR-xQTL_code.R
```

The script contains R code used to reproduce the manuscript figures from processed CSV input files. The input file names are specified in the `read.csv()` calls in the script.

## System requirements

The figure reproduction code was written in R.

Tested environment:

```text
Operating system: Windows
R version: 4.2.1
```

## R package dependencies

The script uses the following R packages:

```r
install.packages(c(
  "circlize",
  "dplyr",
  "forcats",
  "ggbreak",
  "ggplot2",
  "ggpubr",
  "ggridges",
  "ggtern",
  "patchwork",
  "pheatmap",
  "RColorBrewer",
  "survival",
  "survminer",
  "tidyr",
  "VennDiagram"
))
```

The `ComplexHeatmap` and `maftools` packages can be installed from Bioconductor:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("ComplexHeatmap", "maftools"))
```

## Installation

Clone the repository:

```bash
git clone https://github.com/lichenghaoo/ESCC-STR-xQTL.git
cd ESCC-STR-xQTL
```

Install the R package dependencies listed above.

## Data availability

Processed STR-xQTL association results can be queried and downloaded through the ESCC STR-xQTL Atlas:

https://app.su-lab.org/STR-xQTLdb/

Raw WGS, WGBS and RNA-seq data analyzed in this study, from the previously established ESCC multi-omics cohort, have been deposited in the Genome Sequence Archive (GSA) under BioProject accession number PRJCA004501.

Processed ChIP-seq peak data were obtained from ReMap. ESCC GWAS summary statistics were obtained from the GWAS Catalog under accession numbers GCST90271955, GCST90271956 and GCST90018621.

## Usage

Place the required processed CSV files in the working directory, then run:

```bash
Rscript ESCC-STR-xQTL_code.R
```

The script may contain local absolute output paths from the original analysis environment. Before running on a new computer, replace these paths with local output paths, for example:

```r
dir.create("outputs", showWarnings = FALSE)
pdf("outputs/Fig1_panel_B.pdf", width = 6, height = 8)
```

Expected output:

```text
outputs/
|-- Fig1_panel_B.pdf
|-- Fig1_panel_C.pdf
|-- Fig1_panel_D.pdf
`-- additional manuscript figure panel PDFs
```

## Reproducing manuscript figures

The R script is organized by manuscript figure sections. To reproduce a specific figure panel, locate the corresponding section in `ESCC-STR-xQTL_code.R`, ensure that the required processed CSV input files are available, and run the relevant code block.


## License

The figure reproduction code is available for academic review and reuse through the GitHub repository.

## Contact

For questions about code or data, please contact:

```text
lichenghao.bioinformatics@gmail.com
```
