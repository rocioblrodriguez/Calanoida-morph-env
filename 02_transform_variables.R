#!/usr/bin/env Rscript
# =============================================================================
# Script:  02_transform_variables.R
# Author:  Rocio Rodriguez
# Purpose: Standardize environmental predictors and transform morphology.
# Input:   outputs/data_gradients.csv
# Output:  outputs/data_var.csv; outputs/yeojohnson_transformations.rds
# =============================================================================

# -----------------------------------------------------------------------------
# Workflow notes:
# - Starts from the processed organism-level table created in script 01.
# - Standardizes the environmental predictors used later in RDA so they are expressed on comparable scales.
# - Applies a Yeo-Johnson transformation to each selected morphological descriptor while retaining the original raw variables.
# - Saves both the transformed dataset and fitted transformation objects for reproducibility.
# -----------------------------------------------------------------------------


# 1. Settings and packages -----------------------------------------------------
source("R_Scripts/00_setup.R", local = TRUE)
library(readr)
library(dplyr)
library(bestNormalize)

# 2. Read the processed observations ------------------------------------------
dat <- read_csv(file.path(OUT_DIR, "data_gradients.csv"), show_col_types = FALSE)

# 3. Standardize environmental predictors -------------------------------------
# Use all rows in the processed table, matching the original transformation.
# Each standardized variable has mean zero and standard deviation one.
dat <- dat |>
  mutate(across(all_of(ENV_RAW), ~ as.numeric(scale(.x)), .names = "{.col}_std"))

# 4. Transform morphological descriptors --------------------------------------
# Estimate a Yeo–Johnson transformation for each explicitly named descriptor.
# bestNormalize also centers/scales the transformed values. Original columns
# remain available for plots and summaries in their original units.
transformations <- lapply(dat[ALL_MORPH_RAW], yeojohnson, standardize = TRUE)
transformed <- as.data.frame(lapply(transformations, function(x) x$x.t))
names(transformed) <- paste0(ALL_MORPH_RAW, "_YJ")
dat <- bind_cols(dat, transformed)

# 5. Save transformed observations and transformation parameters ---------------
write_csv(dat, file.path(OUT_DIR, "data_var.csv"))
saveRDS(transformations, file.path(OUT_DIR, "yeojohnson_transformations.rds"))
