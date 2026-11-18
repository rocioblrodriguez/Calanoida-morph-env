# =============================================================================
# Script:  helpers.R
# Author:  Rocio Rodriguez
# Purpose: Small shared calculations and figure helpers for the CHP1 dataset.
#          Sources 00_setup.R settings; no input repair or missing-data fallbacks.
# =============================================================================

# -----------------------------------------------------------------------------
# Workflow notes:
# - Contains reusable utility functions shared by multiple analysis and figure scripts.
# - Figure helpers save equivalent output formats and scale organism icons consistently across panels.
# - RDA helpers implement the project-specific abundance transformation, permutation tests, and adjusted-R2 environmental effect calculations.
# - Functions here do not independently read or alter the analysis dataset.
# -----------------------------------------------------------------------------


# 1. Figure exports -----------------------------------------------------------
# Preserve the source plot's dimensions, background, dpi, and device arguments.
# A PNG call creates matching PDF, SVG, and EPS copies. Formats can be restricted for
# a source figure that intentionally uses different PDF and PNG dimensions.
cairo_eps <- function(filename, ...) {
  grDevices::cairo_ps(filename = filename, onefile = FALSE, ...)
}

save_plot_formats <- function(filename, plot, ..., formats = FIGURE_FORMATS) {
  original_args <- list(...)
  stem <- tools::file_path_sans_ext(filename)
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  for (extension in formats) {
    args <- c(list(filename = paste0(stem, ".", extension), plot = plot), original_args)
    if (extension == "pdf") args$device <- grDevices::cairo_pdf
    if (extension == "svg") args$device <- svglite::svglite
    if (extension == "eps") args$device <- cairo_eps
    do.call(ggplot2::ggsave, args)
  }
}

# Map raw ESD to icon width using a common scale across the supplied panels.
icon_sizes <- function(esd, minimum = 0.06, maximum = 0.15) {
  scales::rescale(esd, to = c(minimum, maximum))
}

# 2. RDA helpers --------------------------------------------------------------
# Fit morphology after row multiplication by sqrt(abundance). Only the response
# matrix is multiplied; environmental predictors remain unweighted.
# This preserves the original transformation, followed by rda(scale = TRUE).
fit_morph_rda <- function(dat) {
  Y <- sweep(as.matrix(dat[MORPH_VARIABLES]), 1, sqrt(dat$abundance), FUN = "*")
  vegan::rda(Y ~ ., data = as.data.frame(dat[ENV_STD]), scale = TRUE)
}

# Return full-model and both predictor permutation tests from one fitted RDA.
# All three tests use the same 999 permutations for a reproducible comparison.
# Permutations are unrestricted among rows, as in the supplied analysis.
rda_permutation_tests <- function(fit, n) {
  set.seed(SEED)
  permutations <- permute::shuffleSet(n, nset = N_PERM)
  full <- vegan::anova.cca(fit, permutations = permutations)
  terms <- lapply(c("sequential", "marginal"), function(method) {
    by <- switch(method, sequential = "term", marginal = "margin")
    result <- as.data.frame(vegan::anova.cca(fit, by = by, permutations = permutations))
    result$predictor <- rownames(result)
    result <- result[result$predictor %in% ENV_STD, ]
    data.frame(
      predictor = result$predictor, test = method, Df = result$Df,
      F = result$F, p = result[["Pr(>F)"]], n_permutations = N_PERM
    )
  })
  list(F = full$F[1], p = full[["Pr(>F)"]][1], terms = dplyr::bind_rows(terms))
}

# Return standalone and unique adjusted effects for each predictor.
# Unique effect = adjusted R²(full) - adjusted R²(with that predictor removed).
# Negative estimates are retained in tables; figure links use positive effects.
rda_environment_effects <- function(dat, fit) {
  Y <- sweep(as.matrix(dat[MORPH_VARIABLES]), 1, sqrt(dat$abundance), FUN = "*")
  full_adj <- vegan::RsquareAdj(fit)$adj.r.squared
  dplyr::bind_rows(lapply(ENV_STD, function(variable) {
    reduced <- vegan::rda(
      Y ~ ., data = as.data.frame(dat[setdiff(ENV_STD, variable)]), scale = TRUE
    )
    single <- vegan::rda(Y ~ ., data = as.data.frame(dat[variable]), scale = TRUE)
    reduced_adj <- vegan::RsquareAdj(reduced)$adj.r.squared
    unique_adj <- full_adj - reduced_adj
    data.frame(
      predictor = variable,
      adjR2_single_predictor = vegan::RsquareAdj(single)$adj.r.squared,
      adjR2_unique = unique_adj, reduced_model_adjR2 = reduced_adj,
      full_model_adjR2 = full_adj,
      unique_pct_total = 100 * unique_adj
    )
  }))
}
