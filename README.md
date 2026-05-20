Bayesian differential expression for single-cell RNA-sequencing data via zero-inflated negative bionomial (ZINB) models.

The package wraps Bayesian ZINB modelling for scRNA-seq data stored as Seurat objects using `brms`. It is intended as an alternative to Seurat's default differential expression tests (Wilcoxon rank-sum, DESeq2 or edgeR wrappers).

## Why?
Other than a interesting application of Bayesian modelling in the biologcal field, scRNA-seq count data has two sources of 0 reads:

- **Biological Zeros**: The gene is genuinely not expressed.
- **Technical Dropouts**: Gene expression not captured.

## Installation

SeuratBayesian requires a working Stan installation. You **must** install a working version of `cmdstandr` before executing the workflow. This allows models to compile to C++ on the first run and reduce subsequent computational load.

```r
# install Stan
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
cmdstanr::install_cmdstan()

# install SeuratBayesian
devtools::install_github("aes21/SeuratBayesian")
```

## Quick Start
H

```r
library(Seurat)
library(SeuratData)
library(SeuratBayesian)

# load in example dataset
InstallData("ifnb")
data("ifnb")

#
mono <- subset(ifnb, subset = seurat_annotations == "CD14 Mono")
mono <- NormalizeData(mono)
mono <- FindVariableFeatures(mono, nfeatures = 2000)
gene_of_interest <- head(VariableFeatures(mono), 1) # obvious, but used to display a high magnitude example

VlnPlot_Bayesian(ifnb, feature = gene_of_interest, group.by = "stim", ctr.ident = "CTRL")
```

## Limitations and Additional Comments
As stated, this approach simply aims to
