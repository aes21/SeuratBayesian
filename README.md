# SeuratBayesian

Bayesian differential expression for single-cell RNA-sequencing data via zero-inflated negative binomial (ZINB) models.

The package wraps Bayesian ZINB modelling for scRNA-seq data stored as Seurat objects using `brms`. It is intended for targeted, gene-level characterisation where standard differential expression tools (Wilcoxon rank-sum, DESeq2, edgeR) cannot tell you whether an expression difference reflects a continuous shift across all cells, or a discrete subpopulation switching off entirely.

## Why?
Other than a interesting application of Bayesian modelling in the biological field, standard differential expression tools return a log-fold change and p-value. However, these values summarise an average difference between conditions and cannot decompose the structures behind it.

When a gene shows a large shift in the proportion of cells expressing it - not just in the level of observed expression, a standard test conflates these signals into a single fold-change estimate. The ZINB model separates them explicitly into two components:

- **Negative bionomial**: Expression level in cells where the gene is detected.
- **Zero-inflation**: The probability of a structural zero - a cell where expression is switched off, distinct from a cell with low, but non-zero expression.

The zero-inflation component is fit per condition (`zi ~ condition`). As a result, the model is able to address whether the structural zero probability is changing between conditions.

scRNA-seq count data has two sources of zero reads:

- **Biological Zeros**: The gene is genuinely not expressed.
- **Technical Dropouts**: Gene expression not captured.

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
```

### Standard differential expression approach 
```r
# standard DE approach (Wilcoxon rank sum test)
FindMarkers(mono, feature = "CD14", ident.1 = "STIM", verbose = FALSE)
```

```
             p_val avg_log2FC pct.1 pct.2     p_val_adj
CD14 9.881313e-323  -2.825997 0.253 0.784 1.388621e-318
```

Standard differential expression analysis captures the strong downregulation of CD14 in stimulated monocytes. However, CD14 also exhibits a dramatic drop-rate between the two treatment populations (78% of CTRL cells to just 25% of STIM cells). Here, the fold-change summarises this as an average, but cannot tell if this drop reflects a uniform suppression across the stimulated cells, or a discrete subpopulation of stimulated cells switching CD14 expression off entirely.

### Fit the model to a gene
`sc_fit_bayesian()` returns a standard `brms` fit object, giving direct access to the full suite of brms diagnostics and posterior tools for customisable downstream analysis.

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
```

| Parameter | Description |
|---|---|
| `b_conditionSTIM` | Posterior log fold-change of groups (STIM) vs the control identity (CTRL) in the NB component (e.g., the expression level shift in cells where CD14 is active) |
| `b_zi_Intercept` | Log-odds of a structural zero in CTRL (`plogis(-4.23)` ≈ 1.8%) - CD14 is genuinely expressed |
| `b_zi_conditionSTIM` | Change in log-odds of structural zero in groups (STIM) vs the control identity (CTRL) |
| `shape` | Negative binomial overdispersion; lower values indicate greater count variability |

We can reveal the biological impact of the conditon-level zero-inflation posteriors by converting it to a probability scale:

```r
# convert draws to probability scale
fit_draws <- as_draws_df(fit) |>
  mutate(
    CTRL = plogis(b_zi_Intercept),
    STIM = plogis(b_zi_Intercept + b_zi_conditionSTIM)
  ) |>
  select(CTRL, STIM) |>
  pivot_longer(everything(), names_to = "condition", values_to = "zi_prob")

# summary table
fit_draws |>
  group_by(condition) |>
  summarise(
    median  = median(zi_prob),
    ci_low  = quantile(zi_prob, 0.025),
    ci_high = quantile(zi_prob, 0.975),
    .groups = "drop"
  )
```

```
# A tibble: 2 × 4
  condition median ci_low ci_high
  <chr>      <dbl>  <dbl>   <dbl>
1 CTRL      0.0177 0.0145  0.0214
2 STIM      0.564  0.520   0.611 
```

These results detail how the model interprets control monocytes to show very little "true" zeros (1.8%). As a result, cells in this population exhibited globally lower levels of CD14 expression. Conversely, roughly half of all stimulated monocyte populations are estimated to have structural zero counts (56.4%), suggesting that IFN-β stimulation switches CD14 off entirely in a large subpopulation of cells, while other cells continue to express it, thus driving the original log fold-change observed in standard differential expression analysis.

This distinction matters biologically. Uniform suppression and discrete silencing imply different regulatory mechanisms; only the ZINB model surfaces which is occurring.

### Visualising the results
The whole workflow (model fit to posterior distribution of log fold-change) can be completed using the `VlnPlot_Bayesian()` wrapper function.

```r
VlnPlot_Bayesian(mono, feature = "CD14", group.by = "stim", ctr.ident = "CTRL")
```

![Plot](SeuratBayesian_plot.png)

The top panel details the standard `Seurat::VlnPlot()` output to visual group expression differences across the two treatment groups. The left panel shows the posterior distribution of log fold-change for STIM versus the control reference (fixed at 0). The right panel shows the posterior density of the structural zero probability for each condition. 

Please see the model vignette for justification of the model formula and prior construction.
