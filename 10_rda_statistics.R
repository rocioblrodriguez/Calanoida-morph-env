#!/usr/bin/env Rscript
# =============================================================================
# Script:  10_rda_statistics.R
# Author:  Rocio Rodriguez
# Purpose: Fit RDA models and export overall, sequential, and marginal tests.
# Input:   outputs/Calanoida_data.csv
# Outputs: RDA tables and rda_results.rds, reused by scripts 11 and 12.
# =============================================================================

# -----------------------------------------------------------------------------
# Workflow notes:
# - Fits morphology-environment RDA models from the common clustered Calanoida dataset.
# - Community and regional models are evaluated for Day, Night, and combined DayNight subsets; Region x Cluster node models use daytime observations only.
# - For each model, exports overall fit statistics, constrained-axis variance, sequential and marginal permutation tests, and adjusted-R2 environmental effects.
# - The saved model object is the single source used by scripts 11 and 12, avoiding refitting or inconsistencies between statistics and figures.
# -----------------------------------------------------------------------------


# 1. Settings and packages -----------------------------------------------------
source("R_Scripts/00_setup.R", local = TRUE)
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(vegan)

# 2. Read the common analysis dataset -----------------------------------------
# Regional and community models use Day, Night, and DayNight subsets.
# Cluster × region models use Day only. Groups come from observed data, so an
# absent cluster/region combination does not create an empty model.
dat <- read_csv(file.path(OUT_DIR, "Calanoida_data.csv"), show_col_types = FALSE) |>
  filter(!(Region == "Subarctic" & Cluster == "ST")) |>
  mutate(Region = factor(Region, levels = REGION_LEVELS),
         Cluster = factor(Cluster, levels = CLUSTER_LEVELS))

# 3. Fit and summarize one model ----------------------------------------------
# Return the fitted RDA, observation IDs, model statistics, axis percentages,
# both predictor tests, and environmental effect sizes.
# Sequential tests add Temp, Sal, O2, Fluor in that order. Marginal tests assess
# each predictor after accounting for the other three. Both use 999 permutations.
analyze_rda <- function(subset, diel, approach, group, region = NA_character_, cluster = NA_character_) {
  fit <- fit_morph_rda(subset)
  tests <- rda_permutation_tests(fit, nrow(subset))
  r2 <- RsquareAdj(fit)
  metadata <- tibble(Diel_group = diel, Approach = approach, Group = group,
                     Region = region, Cluster = cluster)
  model <- bind_cols(metadata, tibble(
    n_rows = nrow(subset), full_model_R2 = r2$r.squared,
    full_model_adjR2 = r2$adj.r.squared,
    full_model_pct_total = 100 * r2$r.squared,
    full_model_adj_pct_total = 100 * r2$adj.r.squared,
    F_overall = tests$F, p_overall = tests$p, n_permutations = N_PERM,
    inertia_total = fit$tot.chi, inertia_constrained = fit$CCA$tot.chi,
    inertia_unconstrained = fit$CA$tot.chi
  ))
  predictor_tests <- bind_cols(metadata[rep(1, nrow(tests$terms)), ], tests$terms)
  effects <- rda_environment_effects(subset, fit)
  test_columns <- tests$terms |>
    select(predictor, test, F, p) |>
    pivot_wider(names_from = test, values_from = c(F, p), names_glue = "{.value}_{test}")
  effects <- left_join(effects, test_columns, by = "predictor")
  effects <- bind_cols(metadata[rep(1, nrow(effects)), ], effects)
  eigenvalues <- as.numeric(eigenvals(fit, model = "constrained"))
  axes <- bind_cols(metadata[rep(1, length(eigenvalues)), ], tibble(
    axis = paste0("RDA", seq_along(eigenvalues)), eigenvalue = eigenvalues,
    percent_constrained = 100 * eigenvalues / sum(eigenvalues),
    percent_total = 100 * eigenvalues / fit$tot.chi
  ))
  list(fit = fit, ids = subset$id, model = model, terms = predictor_tests,
       effects = effects, axes = axes)
}

# 4. Run community and regional models ----------------------------------------
results <- list(community = list(), regional = list(), nodes = list())
model_runs <- list()
for (diel in names(DIEL_RUNS)) {
  subset <- filter(dat, D_N %in% DIEL_RUNS[[diel]])
  message("RDA: ", diel, "; community and regional models")
  community <- analyze_rda(subset, diel, "COMMUNITY_WHOLE", "All")
  results$community[[diel]] <- community
  model_runs[[paste(diel, "All", sep = "__")]] <- community
  results$regional[[diel]] <- setNames(lapply(REGION_LEVELS, function(region) {
    analyze_rda(filter(subset, Region == region), diel, "REGION_ENV", region, region)
  }), REGION_LEVELS)
  for (region in REGION_LEVELS) {
    model_runs[[paste(diel, region, sep = "__")]] <- results$regional[[diel]][[region]]
  }
}

# 5. Run daytime cluster × region models --------------------------------------
day_nodes <- dat |>
  filter(D_N == "D") |>
  group_by(Region, Cluster) |>
  group_split()
for (node in day_nodes) {
  region <- as.character(node$Region[1])
  cluster <- as.character(node$Cluster[1])
  node_name <- paste(cluster, region, sep = "_")
  message("RDA: Day; ", node_name)
  result <- analyze_rda(node, "Day", "NODE_ENV", node_name, region, cluster)
  results$nodes[[node_name]] <- result
  model_runs[[paste("Day", node_name, sep = "__")]] <- result
}

# 6. Combine and export tables ------------------------------------------------
# Store signed adjusted effects; zero is not substituted for negative estimates.
# Axis percentages explicitly distinguish constrained from total variance.
collect <- function(component) bind_rows(lapply(model_runs, function(x) x[[component]]))
tables <- list(model_summary = collect("model"), predictor_tests = collect("terms"),
               env_effects = collect("effects"), axis_variance = collect("axes"))
for (name in names(tables)) {
  write_csv(tables[[name]], file.path(OUT_DIR, paste0("RDA_", name, "_all_scales_ALL_DIEL_GROUPS.csv")))
  for (diel in names(DIEL_RUNS)) {
    key <- c(Day = "DAY", Night = "NIGHT", DayNight = "DAY_NIGHT")[[diel]]
    write_csv(filter(tables[[name]], Diel_group == diel),
              file.path(OUT_DIR, paste0("RDA_", name, "_all_scales_", key, ".csv")))
  }
}

# Daytime node tables provide a compact supplementary-table source and the
# exact values used by the environmental-link figure.
write_csv(filter(tables$model_summary, Approach == "NODE_ENV"),
          file.path(OUT_DIR, "RDA_cluster_region_DAY_model_summary.csv"))
write_csv(filter(tables$predictor_tests, Approach == "NODE_ENV"),
          file.path(OUT_DIR, "RDA_cluster_region_DAY_predictor_tests_both.csv"))
write_csv(filter(tables$env_effects, Approach == "NODE_ENV"),
          file.path(OUT_DIR, "RDA_cluster_region_DAY_environment_effects.csv"))
saveRDS(results, file.path(OUT_DIR, "rda_results.rds"))
