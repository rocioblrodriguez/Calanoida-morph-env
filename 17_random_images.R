#!/usr/bin/env Rscript
# =============================================================================
# Script:  17_random_images.R
# Author:  Rocio Rodriguez
# Purpose: Select 20 random organism images per Region × Cluster combination.
# Input:   outputs/Calanoida_data.csv and the EcoTaxa JPG directory in setup.
# Outputs: Copied JPGs, a selection manifest, and group counts in outputs/.
# =============================================================================

# 1. Settings and packages -----------------------------------------------------
source("R_Scripts/00_setup.R", local = TRUE)
library(readr)
library(dplyr)
library(tidyr)

IMAGES_PER_GROUP <- 20L
RANDOM_IMAGE_SEED <- SEED
SELECTION_DIR <- file.path(OUT_DIR, "random_images_20_per_cluster_region")
MANIFEST_FILE <- file.path(OUT_DIR, "random_images_20_per_cluster_region.csv")
COUNTS_FILE <- file.path(OUT_DIR, "random_images_20_per_cluster_region_counts.csv")

# 2. Match clustered observations to available JPGs ----------------------------
dat <- read_csv(file.path(OUT_DIR, "Calanoida_data.csv"), show_col_types = FALSE) |>
  filter(
    !is.na(Region), !is.na(Cluster), !is.na(id),
    !(Region == "Subarctic" & Cluster == "ST")
  ) |>
  distinct(id, .keep_all = TRUE) |>
  mutate(
    Region = factor(Region, levels = REGION_LEVELS),
    Cluster = factor(Cluster, levels = CLUSTER_LEVELS),
    image_path = file.path(IMAGES_DIR, paste0(id, ".jpg"))
  )

available <- dat[file.exists(dat$image_path), , drop = FALSE]
if (!nrow(available)) {
  stop("No organism JPGs matched the id values. Check IMAGES_DIR in 00_setup.R: ",
       IMAGES_DIR)
}

missing_images <- nrow(dat) - nrow(available)
if (missing_images > 0L) {
  warning(missing_images, " clustered observations had no matching JPG and were excluded.")
}

# 3. Draw an equal-probability sample without replacement ----------------------
set.seed(RANDOM_IMAGE_SEED)
selected <- available |>
  group_by(Region, Cluster) |>
  group_modify(~ slice_sample(.x, n = min(IMAGES_PER_GROUP, nrow(.x)))) |>
  ungroup() |>
  arrange(Region, Cluster, id)

# 4. Copy selected images into Region/Cluster folders --------------------------
# This directory contains generated copies only, so recreate it on each run to
# prevent images from an older selection remaining beside the current sample.
if (dir.exists(SELECTION_DIR)) unlink(SELECTION_DIR, recursive = TRUE)
dir.create(SELECTION_DIR, recursive = TRUE, showWarnings = FALSE)

selected <- selected |>
  mutate(
    destination_dir = file.path(SELECTION_DIR, as.character(Region), as.character(Cluster)),
    destination_path = file.path(destination_dir, paste0(id, ".jpg"))
  )
invisible(lapply(unique(selected$destination_dir), dir.create,
                 recursive = TRUE, showWarnings = FALSE))

copied <- file.copy(selected$image_path, selected$destination_path, overwrite = TRUE)
if (!all(copied)) {
  stop("Failed to copy ", sum(!copied), " selected organism JPGs.")
}

# 5. Save the selected observations and counts ---------------------------------
write_csv(
  selected |>
    select(Region, Cluster, id, abundance, Cruise, Station, D_N,
           image_path, destination_path),
  MANIFEST_FILE
)

available_counts <- available |>
  count(Region, Cluster, name = "available_images")
selected_counts <- selected |>
  count(Region, Cluster, name = "selected_images")

count_summary <- expand_grid(
  Region = factor(REGION_LEVELS, levels = REGION_LEVELS),
  Cluster = factor(CLUSTER_LEVELS, levels = CLUSTER_LEVELS)
) |>
  left_join(available_counts, by = c("Region", "Cluster")) |>
  left_join(selected_counts, by = c("Region", "Cluster")) |>
  mutate(
    requested_images = IMAGES_PER_GROUP,
    available_images = replace_na(available_images, 0L),
    selected_images = replace_na(selected_images, 0L)
  ) |>
  arrange(Region, Cluster)
write_csv(count_summary, COUNTS_FILE)

message(
  "Random image selection complete: ", nrow(selected), " JPGs copied to ",
  SELECTION_DIR, ". Seed: ", RANDOM_IMAGE_SEED, "."
)
