#!/usr/bin/env Rscript
# =============================================================================
# Script:  04_morphology_correlations.R
# Author:  Rocio Rodriguez
# Purpose: Calculate correlations among all transformed ZooProcess descriptors.
# Input:   outputs/data_var_clusters.csv
# Output:  outputs/cor_all_morph_vars_YJ.csv
# =============================================================================

# 1. Settings and packages -----------------------------------------------------
source("R_Scripts/00_setup.R", local = TRUE)
library(readr)
library(dplyr)
library(tibble)

# 2. Select validated NH1208 Calanoida -----------------------------------------
dat <- read_csv(file.path(OUT_DIR, "data_var_clusters.csv"), show_col_types = FALSE) |>
  filter(
    Cruise == CRUISE,
    annotation_category == CATEGORY,
    annotation_status == STATUS,
    !(Region == "Subarctic" & Cluster == "ST")
  )

# 3. Calculate Pearson correlations -------------------------------------------
# This diagnostic includes the full descriptor set, including position and
# bounding-box variables that are excluded from the 43-variable main analysis.
correlations <- cor(dat[paste0(ALL_MORPH_RAW, "_YJ")], method = "pearson")

# 4. Save one row per descriptor ----------------------------------------------
write_csv(
  rownames_to_column(as.data.frame(correlations), "variable"),
  file.path(OUT_DIR, "cor_all_morph_vars_YJ.csv")
)
