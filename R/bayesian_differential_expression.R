#' Fit ZINB Bayesian modelling for Seurat gene expression.
#'
#' @param object A Seurat object.
#' @param feature Genes to test.
#' @param group.by Regroup cells into a different identity class prior to
#' performing differential expression (see example); "ident" to use Idents
#' @param ctr.ident The control identity group.
#' @param iter A positive integer specifying the number of iterations for each
#' chain (including warmup). The default is 1000.
#' @param warmup 	A positive integer specifying the number of warmup (aka
#' burnin) iterations per chain.
#' @param chains A positive integer specifying the number of Markov chains.
#' @param cores Number of cores to use when executing the chains in parallel,
#' defaults to \code{parallel::detectCores}. See: \code{\link[rstan]{sampling}}.
#' @param seed The seed for random number generation.
#' @param ... Additional arguments passed to \code{\link[brms]{brm}}.
#'
#' @seealso \code{\link[brms]{brm}} for the full list of sampling arguments.
#'
#' @import brms
#' @import dplyr
#' @import Seurat
#' @importFrom Matrix colSums
#' @importFrom parallel detectCores
#'
#' @return data.frame with ZINB model summary metrics.
#' 
#' @examples
#' \dontrun{
#' data("pbmc_small")
#' fit <- sc_fit_bayesian(object = pbmc_small, feature = "LYZ", group.by = "groups", ctr.ident = "g1")
#' }
#'
#' @export
#'
sc_fit_bayesian <- function(
  object,
  feature,
  group.by,
  ctr.ident,
  iter = 1000L,
  warmup = 500L,
  chains = 4L,
  cores = parallel::detectCores(logical = FALSE),
  seed = 42L,
  ...
) {
  # validate parameters
  stopifnot(
    inherits(object, "Seurat"),
    feature %in% rownames(object),
    length(feature) == 1L,
    group.by %in% colnames(object@meta.data),
    ctr.ident %in% object@meta.data[[group.by]]
  )

  # ensure formatting
  count_mtx <- Seurat::GetAssayData(
    object,
    assay = Seurat::DefaultAssay(object), 
    layer = "counts"
  )

  # empty count matrix due to incorrect Seurat layer
  if (ncol(count_mtx) == 0 || nrow(count_mtx) == 0) {
    stop("The count matrix is empty, you must ensure the 'DefaultAssay' of this Seurat object has a 'counts' layer.")
  }

  # extract data from Seurat object
  lib_size <- Matrix::colSums(count_mtx)
  gene_data <- data.frame(
    cell_id = colnames(object),
    counts = as.numeric(count_mtx[feature, ]),
    condition = object@meta.data[[group.by]],
    log_lib_size = log(lib_size)
  )

  # set control group as first identity
  gene_data$condition <- relevel(factor(gene_data$condition), ref = ctr.ident)

  # fit zinb model
  zinb_formula <- brms::bf(counts ~ condition +
                             offset(log_lib_size), zi ~ condition)
  fit <- brms::brm(
    formula = zinb_formula,
    data = gene_data,
    family = brms::zero_inflated_negbinomial(),
    algorithm = "meanfield",
    prior = c(
      brms::prior(normal(0, 2), class = "b"),
      brms::prior(normal(0, 1.5), class = "Intercept", dpar = "zi")
    ),
    backend = "cmdstanr",
    stan_model_args = list(stanc_options = list("O1")),
    ...
  )

  return(fit)
}

#' Plot Bayesian fit predictions.
#'
#' @param object A Seurat object.
#' @param feature Genes to test.
#' @param group.by Regroup cells into a different identity class prior to
#' performing differential expression (see example); "ident" to use Idents
#' @param ctr.ident The control identity group.
#' @param cols Colors to use for plotting.
#'
#' @import dplyr
#' @import ggplot2
#' @import tidyr
#' @importFrom posterior as_draws_df
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' data("pbmc_small")
#' VlnPlot_Bayesian(object = pbmc_small, feature = "LYZ", group.by = "groups", ctr.ident = "g1")
#' }
#'
#' @export
#'
VlnPlot_Bayesian <- function(
  object,
  feature,
  group.by,
  ctr.ident,
  cols = NULL
) {
  # run fit
  fit <- sc_fit_bayesian(
    object = object,
    feature = feature,
    group.by = group.by,
    ctr.ident = ctr.ident
  )

  # generate formatted draws object
  draws <- as_draws_df(fit) |>
    select(starts_with("b_condition")) |>
    pivot_longer(cols = everything(),
                 names_to = "condition",
                 values_to = "lfc") |>
    mutate(condition = gsub("b_condition", "", condition))

  # add control reference (LFC = 0 by definition)
  control_draws <- tidyr::tibble(condition = ctr.ident, lfc = 0)
  draws <- dplyr::bind_rows(control_draws, draws)

  # compute credible intervals for annotation
  ci_df <- draws |>
    group_by(condition) |>
    summarise(
      median = median(lfc),
      ci_low  = quantile(lfc, 0.025),
      ci_high = quantile(lfc, 0.975),
      .groups = "drop"
    )

  ggplot(draws, aes(x = condition, y = lfc, fill = condition)) +
    geom_violin(trim = TRUE,
                scale = "width",
                alpha = 0.8) +
    geom_hline(yintercept = 0,
               linetype = "dashed",
               colour = "grey40") +
    geom_pointrange(
      data = ci_df,
      aes(y = median, ymin = ci_low, ymax = ci_high),
      colour = "transparent",
      size = 0.5,
      linewidth = 0.8
    ) +
    labs(
      title = feature,
      x = NULL,
      y = "Log Fold-Change"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}
