#!/usr/bin/env Rscript
# =============================================================================
# Script:  03_cluster_morphology.R
# Author:  Rocio Rodriguez
# Purpose: Fit the six-morphotype solution and the two-cluster comparison.
# Input:   outputs/data_var.csv
# Outputs: Cluster assignments, model objects, and regional summaries.
# =============================================================================

# 1. Settings and packages -----------------------------------------------------
source("R_Scripts/00_setup.R", local = TRUE)
library(readr)
library(dplyr)

# 2. Select observations -------------------------------------------------------
# Fit the clusters using all validated NH1208 Calanoida. Subarctic ST can only
# be identified after cluster labels are assigned, so it is removed below before
# any clustered dataset, summary, PCA, behavior, RDA, or figure is saved.
dat <- read_csv(file.path(OUT_DIR, "data_var.csv"), show_col_types = FALSE) |>
  filter(Cruise == CRUISE, annotation_category == CATEGORY, annotation_status == STATUS)

# 3. Fit six clusters using the final 20-variable set --------------------------
# YJ variables were standardized in script 02. As in the supplied six-cluster
# analysis, use those values directly and give each image equal influence.
set.seed(SEED)
km_six <- kmeans(dat[CLUSTERING_VARIABLES], centers = 6, nstart = KMEANS_NSTART, iter.max = 300)
dat$Cluster_id <- km_six$cluster

# 4. Assign descriptive labels ------------------------------------------------
# Order cluster centroids by mean grayscale intensity: the upper three are
# transparent. Within transparent and dark groups, order mean ESD from small
# to large. These labels summarize clusters; they do not define their fit.
cluster_stats <- dat |>
  group_by(Cluster_id) |>
  summarise(mean_transparency = mean(mean_YJ), mean_size = mean(esd_YJ), .groups = "drop") |>
  arrange(desc(mean_transparency)) |>
  mutate(transparency_cat = if_else(row_number() <= 3, "transparent", "dark")) |>
  group_by(transparency_cat) |>
  arrange(mean_size, .by_group = TRUE) |>
  mutate(cluster_label = paste0(c("S", "M", "L")[row_number()],
                                if_else(transparency_cat == "transparent", "T", "D"))) |>
  ungroup()
label_map <- setNames(cluster_stats$cluster_label, cluster_stats$Cluster_id)
dat$Cluster <- factor(unname(label_map[as.character(dat$Cluster_id)]), levels = CLUSTER_LEVELS)

# 5. Fit the independent two-cluster comparison -------------------------------
# Preserve the original comparison: k = 2 uses all 43 YJ descriptors,
# standardized within this subset, without abundance multiplication.
# S and L are assigned by the unweighted mean raw ESD of the two groups.
two_cluster_matrix <- scale(dat[MORPH_VARIABLES])
set.seed(SEED)
km_two <- kmeans(two_cluster_matrix, centers = 2, nstart = KMEANS_NSTART, iter.max = 1000)
size_by_group <- tapply(dat$esd, km_two$cluster, mean)
two_labels <- setNames(c("S", "L"), names(sort(size_by_group)))
dat$Cluster2_id <- km_two$cluster
dat$Cluster2 <- unname(two_labels[as.character(km_two$cluster)])

# Remove ST only within the Subarctic region from every downstream analysis.
exclude_subarctic_st <- !is.na(dat$Region) & !is.na(dat$Cluster) &
  as.character(dat$Region) == "Subarctic" &
  as.character(dat$Cluster) == "ST"
message("Excluding ", sum(exclude_subarctic_st, na.rm = TRUE),
        " Subarctic ST observations from all downstream analyses and plots.")
dat <- dat[!exclude_subarctic_st, , drop = FALSE]

# 6. Summarize regional representation ----------------------------------------
cluster_region_counts <- dat |>
  group_by(Region, Cluster) |>
  summarise(
    n_rows = n(),
    n_D = sum(D_N == "D", na.rm = TRUE),
    n_N = sum(D_N == "N", na.rm = TRUE),
    n_D_minus_N = n_D - n_N,
    sum_abundance = sum(abundance, na.rm = TRUE),
    abundance_D = sum(abundance[D_N == "D"], na.rm = TRUE),
    abundance_N = sum(abundance[D_N == "N"], na.rm = TRUE),
    abundance_D_minus_N = abundance_D - abundance_N,
    mean_abundance = mean(abundance, na.rm = TRUE),
    .groups = "drop"
  )
selected_variables <- data.frame(
  variable = CLUSTERING_VARIABLES,
  category = rep(names(CLUSTERING_GROUPS), lengths(CLUSTERING_GROUPS))
)

# 7. Save assignments, settings, and fitted models -----------------------------
write_csv(dat, file.path(OUT_DIR, "data_var_clusters.csv"))
write_csv(cluster_stats, file.path(OUT_DIR, "cluster_stats.csv"))
cluster_region_file <- file.path(OUT_DIR, "cluster_rows_per_region.csv")
write_csv(cluster_region_counts, cluster_region_file)

# Confirm that this run replaced any older five-column version of the table.
expected_region_columns <- c(
  "Region", "Cluster", "n_rows", "n_D", "n_N", "n_D_minus_N",
  "sum_abundance", "abundance_D", "abundance_N",
  "abundance_D_minus_N", "mean_abundance"
)
written_region_columns <- names(readr::read_csv(cluster_region_file, n_max = 0, show_col_types = FALSE))
if (!identical(written_region_columns, expected_region_columns)) {
  stop("cluster_rows_per_region.csv was not written with the expected D/N columns.")
}
message("Updated cluster_rows_per_region.csv with D/N observation and abundance differences.")
write_csv(selected_variables, file.path(OUT_DIR, "selected_20_clustering_variables.csv"))
saveRDS(
  list(six = km_six, two = km_two, six_labels = label_map, two_labels = two_labels,
       six_variables = CLUSTERING_VARIABLES, two_variables = MORPH_VARIABLES,
       two_center = attr(two_cluster_matrix, "scaled:center"),
       two_scale = attr(two_cluster_matrix, "scaled:scale"), seed = SEED),
  file.path(OUT_DIR, "kmeans_models.rds")
)
