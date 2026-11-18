#!/usr/bin/env Rscript
# =============================================================================
# Script:  00_run_all.R
# Author:  Rocio Rodriguez
# Purpose: Reproduce the CHP1 tables and figures in their dependency order.
# Input:   Raw files and archived diagnostics listed in README.md.
# Outputs: outputs/ and plots/, plus a record of R and package versions.
# =============================================================================

# 1. Execution order ----------------------------------------------------------
# The scripts are numbered in dependency order. This explicit list excludes
# this runner, setup, and optional installation steps.
scripts <- c(
  "01_prepare_data.R", "02_transform_variables.R", "03_cluster_morphology.R",
  "04_morphology_correlations.R", "05_pca_statistics.R", "06_pca_plots.R",
  "07_pca_two_clusters.R", "08_size_transparency.R", "09_behavior.R",
  "10_rda_statistics.R", "11_rda_images.R", "12_centroid_environment_links.R",
  "13_butterfly_profiles.R", "14_silhouette_maxima.R", "15_sampling_map.R",
  "16_environmental_sections.R", "17_random_images.R"
)

# 2. Run the workflow ---------------------------------------------------------
# A separate environment for each script keeps its temporary objects local.
# An analysis error stops execution directly instead of producing partial runs.
for (script in scripts) {
  message("Running ", script)
  source(file.path("R_Scripts", script), local = new.env(parent = globalenv()))
}

# 3. Record the computational environment -------------------------------------
capture.output(sessionInfo(), file = file.path("outputs", "sessionInfo.txt"))
write.csv(
  as.data.frame(installed.packages()[, c("Package", "Version")]),
  file.path("outputs", "package_versions.csv"), row.names = FALSE
)
message("CHP1 workflow complete.")
