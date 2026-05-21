# SeuratBayesian

Bayesian differential expression for single-cell RNA-sequencing data via zero-inflated negative binomial (ZINB) models.

The package wraps Bayesian ZINB modelling for scRNA-seq data stored as Seurat objects using `brms`. It is intended for targeted, gene-level characterisation where standard differential expression screening (Wilcoxon rank-sum, DESeq2, edgeR) returns ambiguous results or where understanding the mechanism behind a significant expression difference is important.

## Why?
Other than a interesting application of Bayesian modelling in the biological field, the standard differential expression tools return a log-fold change and p-value. However, these values do not tell you why the counts differ between conditions. scRNA-seq count data has two sources of 0 reads:

- **Biological Zeros**: The gene is genuinely not expressed.
- **Technical Dropouts**: Gene expression not captured.

Where a gene shows a large shift in the proportion of cells expressing it - not just in the level of observed expression, standard differential expression analysis conflates these two signals. The proposed ZINB model separates them with a negative binomial component for expression level, and a zero-inflation component for the probability of a structural zero. This decomposition gives a more complete picture of what is actually changing between conditions.

## Installation

SeuratBayesian requires a working Stan installation. You **must** install `cmdstanr` before executing the workflow. This allows models to compile to C++ on the first run and reduce subsequent computational load.

```r
# install Stan
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
cmdstanr::install_cmdstan()

# install SeuratBayesian
devtools::install_github("aes21/SeuratBayesian")
```

## Quick Start
The following example uses the `ifnb` dataset from the SeuratData package. Using a subset of CD14+ monocytes, we examine what the ZINB model reveals about the expression of CD14 under stimulation treatment.

```r
library(Seurat)
library(SeuratBayesian)
library(SeuratData)

InstallData("ifnb")
data("ifnb")

# subset for CD14+ monocyte cell group
mono <- subset(ifnb, subset = seurat_annotations == "CD14 Mono")

# standard DE approach (Wilcoxon rank sum test)
FindMarkers(mono, feature = "CD14", ident.1 = "STIM", verbose = FALSE)
```

```
             p_val avg_log2FC pct.1 pct.2     p_val_adj
CD14 9.881313e-323  -2.825997 0.253 0.784 1.388621e-318
```

In this scenario, standard differential analysis captures the direction and size of CD14 expression changes (near-zero p-value and large log fold-change across the stimulated group), but cannot tell you whether the significant drop in detection rates (78% in CTRL to 25% in STIM) reflects suppressed transcription across all cells, a subpopulation switching off, or increased technical dropout as a result of the treatment. The posterior characterises both the expression level shift and the zero-inflation structure simultaneously, providing a more informative answer than a p-value alone.

### Fit the model to a gene
The `sc_fit_bayesian()` function returns a standard `brms` fit object, giving direct access to the full suite of brms diagnostics and posterior tools for custom downstream analysis.

```r
# fit model
fit <- sc_fit_bayesian(
  object = mono,
  feature = "CD14",
  group.by = "stim",
  ctr.ident = "CTRL"
)

# downstream brms tools to evaluate fit
brms::posterior_summary(fit)
```

```
                   Estimate  Est.Error        Q2.5      Q97.5
b_Intercept      -6.4641956 0.04050887  -6.5386461 -6.3821472
b_zi_Intercept   -4.2316024 0.65323779  -5.4748927 -2.9132264
b_conditionSTIM  -1.9981520 0.06535874  -2.1255479 -1.8734239
shape             0.7940529 0.03202204   0.7357693  0.8598465
Intercept        -7.4476968 0.02589520  -7.4989307 -7.3960797
Intercept_zi     -4.2316024 0.65323779  -5.4748927 -2.9132264
lprior          -11.0757485 1.22774904 -13.7458231 -8.8887714
lp_approx__      -1.9977503 1.40918513  -5.5777273 -0.2784322
```

- **b_Intercept**: Baseline log expression level in the control identity, in this case reflecting the sparisty of CD14 counts.
- **b_zi_Intercept**: The probability of a structural 0 from the ZINB component across the conditions. A 1.4% low value (`plogis(-4.23)`) indicates the 0 reads of CD14 counts are largely explained by a genuine reduction of expression, as opposed to dropouts and reflects a true, biological suppression.
- **b_conditionSTIM**: Posterior log fold-change of STIM vs. CTRL, mimics the observed downregulation observed using the Wilcoxon rank sum test.

### Visualising posterior distributions
The whole workflow (model fit to posterior distribution of log fold-change) can be completed using the `VlnPlot_Bayesian()` wrapper function.

```r
VlnPlot_Bayesian(mono, feature = gene_of_interest, group.by = "stim", ctr.ident = "CTRL")
```

![Plot](SeuratBayesian_plot.png)

The violin plot shows the posterior distribution of log fold-change for each group against the defined control. The control group is shown with a log fold-change fixed at 0 (as the reference level). Please see the model vignette for justification of the model formula and prior construction.
