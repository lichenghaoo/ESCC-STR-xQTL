# ESCC-STR-xQTL

Figure reproduction code for the multi-omics STR-xQTL study of esophageal squamous cell carcinoma.

This repository accompanies the manuscript:

> Multi-omics analysis of short tandem repeat variation identifies candidate causal risk loci in esophageal squamous cell carcinoma

## Scope

This repository contains custom R scripts for reproducing manuscript figures from processed analysis outputs. It is intended to support code review and figure reproduction for the manuscript.

This repository is not a standalone software package and is not a complete end-to-end pipeline for raw WGS, WGBS or RNA-seq processing from FASTQ files.

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

No non-standard hardware is required for running the figure reproduction scripts on processed data. A standard desktop or laptop computer is sufficient for the plotting scripts.

Full raw sequencing data processing is outside the scope of this repository and may require a high-performance computing environment.

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

The code uses processed outputs generated from the multi-omics STR-xQTL analyses described in the manuscript. It does not rerun STR genotyping, raw sequencing alignment, molecular phenotype quantification, STR-xQTL mapping, fine-mapping or colocalization from raw sequencing files.

## Software and analysis tools used in the manuscript

The full manuscript analyses used established public tools, including ExpansionHunter v5.0.0, GangSTR v2.5.0, TRTools v6.0.2, EnsembleTR v1.0.0, STAR v2.7.11, RSEM v1.3.1, LeafCutter v0.2.9, DaPars v2.1, Bismark v0.10.1, LDAK v6.1, ChIPseeker v1.24.0, bedtools v2.31.1, bmediatR v0.1.3, susieR v0.12.35, coloc v5.2.3, IGV v2.18.0 and Cytoscape v3.10.4.

The present repository focuses on custom R scripts for figure reproduction from processed analysis outputs.

## License

The figure reproduction code is available for academic review and reuse through the GitHub repository.

## Contact

For questions about code or data, please contact:

```text
lichenghao.bioinformatics@gmail.com
```
