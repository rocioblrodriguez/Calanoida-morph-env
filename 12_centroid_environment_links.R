#!/usr/bin/env Rscript
# =============================================================================
# Script:  12_centroid_environment_links.R
# Author:  Rocio Rodriguez
# Purpose: Compare morphological centroids and show daytime RDA effects.
# Inputs:  Calanoida_data.csv and rda_results.rds from scripts 05 and 10.
# Outputs: Correlation figure, matrices, both predictor tests, and link tables.
# =============================================================================

# -----------------------------------------------------------------------------
# Workflow notes:
# - Calculates abundance-weighted Region x Cluster morphological centroids from the 43 transformed descriptors.
# - Uses centroid correlations to describe morphological similarity among occupied nodes.
# - Reuses the daytime node-level RDA statistics from script 10 to connect environmental effects with those morphological nodes.
# - Exports the correlation matrices, predictor-test tables, environmental-link data, and the final combined visualization.
# -----------------------------------------------------------------------------


# 1. Settings and packages -----------------------------------------------------
source("R_Scripts/00_setup.R", local = TRUE)
PLOT_DIR <- PLOT_DIRS[["Centroid_environment_links"]]
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(Hmisc)
library(ggcor)
library(ggplot2)
ENV_LINK_TEST <- "marginal" # Both tests are exported; this one selects plotted links.
ENV_LINK_ALPHA <- 0.05
TOP_ENV_LINKS <- 2
LINK_WIDTH_RANGE <- c(0.3, 3.5)
LINK_LEGEND_BREAKS <- c(0, 5, 10)

# 2. Calculate abundance-weighted morphological centroids -----------------------
# Both day and night contribute to each region × cluster centroid. Correlations
# compare the 43-descriptor centroid profiles; the RDA links use daytime data.
dat <- read_csv(file.path(OUT_DIR, "Calanoida_data.csv"), show_col_types = FALSE) |>
  filter(
    Cruise == CRUISE,
    annotation_category == CATEGORY,
    annotation_status == STATUS,
    Region %in% REGION_LEVELS,
    Cluster %in% NODE_CLUSTER_LEVELS,
    !(Region == "Subarctic" & Cluster == "ST"),
    is.finite(abundance),
    abundance > 0
  ) |>
  mutate(
    Region = factor(Region, levels = REGION_LEVELS),
    Cluster = factor(trimws(as.character(Cluster)), levels = NODE_CLUSTER_LEVELS)
  ) |>
  drop_na(abundance, D_N, all_of(ENV_RAW), all_of(MORPH_VARIABLES))

if (!nrow(dat)) {
  stop("No complete positive-abundance rows are available for the centroid plot.")
}

centroids <- dat |>
  group_by(Region, Cluster) |>
  summarise(
    across(all_of(MORPH_VARIABLES), ~ weighted.mean(.x, abundance, na.rm = TRUE)),
    total_abundance = sum(abundance, na.rm = TRUE),
    n_objects = n(),
    .groups = "drop"
  ) |>
  mutate(node = paste(Cluster, Region, sep = "_")) |>
  arrange(Region, Cluster)
node_matrix <- centroids |>
  select(node, all_of(MORPH_VARIABLES)) |>
  column_to_rownames("node") |>
  as.matrix() |>
  t()

# Preserve the original node arrangement while excluding ST only in Subarctic.
# Transition and Subtropical retain ST and all other observed clusters.
expected_nodes <- unlist(lapply(
  REGION_LEVELS,
  function(region) paste(NODE_CLUSTER_LEVELS, region, sep = "_")
))
expected_nodes <- setdiff(expected_nodes, "ST_Subarctic")
missing_nodes <- setdiff(expected_nodes, colnames(node_matrix))
if (length(missing_nodes)) {
  stop(
    "The expected 17-node layout cannot be produced because these nodes are missing: ",
    paste(missing_nodes, collapse = ", ")
  )
}
node_matrix <- node_matrix[, expected_nodes, drop = FALSE]
write_csv(centroids, file.path(OUT_DIR, "morphological_centroids_cluster_region.csv"))

# 3. Export both predictor tests from the same daytime RDA models ---------------
# Script 10 fitted and tested these exact models. Reuse its results to prevent
# p-values or effect definitions from drifting between the table and figure.
rda <- readRDS(file.path(OUT_DIR, "rda_results.rds"))
tests <- bind_rows(lapply(rda$nodes, function(x) x$terms))
effects <- bind_rows(lapply(rda$nodes, function(x) x$effects))
write_csv(tests, file.path(OUT_DIR, "node_RDA_predictor_tests_sequential_and_marginal.csv"))
links <- effects |>
  mutate(node = Group, environment = sub("_std$", "", predictor),
         p.value = .data[[paste0("p_", ENV_LINK_TEST)]],
         effect_pct = unique_pct_total, test_used_for_plot = ENV_LINK_TEST)
write_csv(links, file.path(OUT_DIR, "node_environment_links_all_unique_effects_abundance_weighted.csv"))

# Rank significant positive effects independently within each region × cluster.
# The signed unique adjusted R² values remain available in the complete table.
selected_links <- links |>
  filter(p.value <= ENV_LINK_ALPHA, effect_pct > 0) |>
  group_by(Region, Cluster) |>
  arrange(desc(effect_pct), p.value, environment, .by_group = TRUE) |>
  slice_head(n = TOP_ENV_LINKS) |>
  mutate(link_rank = row_number()) |>
  ungroup() |>
  mutate(.row = unname(c(Temp = "Temperature", Sal = "Salinity", O2 = "Oxygen", Fluor = "Fluorescence")[environment]),
         .col = node)
write_csv(selected_links, file.path(OUT_DIR, "node_environment_links_top2_significant_abundance_weighted.csv"))

# 4. Calculate Pearson correlations and their reported p-values ----------------
correlations <- rcorr(node_matrix, type = "pearson")
R <- correlations$r
P <- correlations$P
write.csv(R, file.path(OUT_DIR, "node_correlations_r_full_abundance_weighted.csv"))
write.csv(P, file.path(OUT_DIR, "node_correlations_p_full_abundance_weighted.csv"))
R_significant <- R
R_significant[P >= 0.05] <- NA_real_
diag(R_significant) <- 1
write.csv(R_significant, file.path(OUT_DIR, "node_correlations_significant_matrix_abundance_weighted.csv"))

# The correlation significance tests use descriptors as observations, as in the
# supplied analysis. The p-values are unadjusted; interpretation should account
# for correlations among descriptors and the number of pairwise comparisons.
nodes <- colnames(node_matrix)
N <- length(nodes)
pairs <- as.data.frame(as.table(R)) |>
  rename(node1 = Var1, node2 = Var2, r = Freq) |>
  left_join(as.data.frame(as.table(P)) |> rename(node1 = Var1, node2 = Var2, p = Freq),
            by = c("node1", "node2")) |>
  mutate(row_index = match(node1, nodes), col_index = match(node2, nodes),
         sig = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")) |>
  filter(col_index >= row_index)
write_csv(filter(pairs, sig != ""), file.path(OUT_DIR, "node_correlations_significant_long_abundance_weighted.csv"))

# Original graphical inputs, using the exact signed effects exported above.
vars <- nodes
cluster_levels <- NODE_CLUSTER_LEVELS
region_levels <- REGION_LEVELS
env_links_top <- selected_links
link_size_upper <- max(c(LINK_LEGEND_BREAKS, selected_links$effect_pct))
env_palette_full <- c(
  "Temperature (°C)" = "#E64B35", "Oxygen (μmol kg⁻¹)" = "#4DBBD5",
  "Salinity (PSU)" = "#fc8d62", "Fluorescence (mg m⁻³)" = "#00A087"
)
# ---- 6) P-STARS FOR PLOTTING ------------------------------------------------

ix <- setNames(seq_along(vars), vars)

p_star_df <- as.data.frame(as.table(P), stringsAsFactors = FALSE) %>%
  rename(row_lab = Var1, col_lab = Var2, p = Freq) %>%
  mutate(
    row_i = ix[row_lab],
    col_i = ix[col_lab],
    keep  = col_i >= row_i,
    x     = col_i,
    y     = N - row_i + 1,
    sig   = case_when(
      is.na(p)  ~ NA_character_,
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE      ~ NA_character_
    )
  ) %>%
  filter(keep, !is.na(sig))

# ---- 7) LABELS + REGION LINES -----------------------------------------------

node_info <- tibble(node = vars) %>%
  tidyr::separate(
    node,
    into = c("Cluster", "Region"),
    sep = "_",
    remove = FALSE
  ) %>%
  mutate(
    Cluster = factor(Cluster, levels = cluster_levels),
    Region  = factor(Region,  levels = region_levels),
    x       = seq_along(node),
    y       = N - x + 1L
  )

region_blocks_x <- node_info %>%
  group_by(Region) %>%
  summarise(xmid = mean(range(x)), .groups = "drop")

region_blocks_y <- node_info %>%
  group_by(Region) %>%
  summarise(ymid = mean(range(y)), .groups = "drop")

clabs <- tibble(
  x   = seq_along(vars),
  lab = sub("_.*", "", vars)
)

ylabs <- tibble(
  y   = N - seq_along(vars) + 1L,
  lab = sub("_.*", "", vars)
)

top_y   <- N + 0.9
right_x <- N + 0.9

# Region delineation lines on upper triangle only.
# Boundaries follow the retained nodes after excluding Subarctic ST.
region_breaks <- node_info %>%
  group_by(Region) %>%
  summarise(
    xmin = min(x),
    xmax = max(x),
    ymin = min(y),
    ymax = max(y),
    .groups = "drop"
  ) %>%
  arrange(Region)

break_after <- region_breaks$xmax[-nrow(region_breaks)] + 0.5

vseg <- tibble(
  x    = break_after,
  xend = break_after,
  y    = N - break_after + 1,
  yend = N + 0.5
)

hseg <- tibble(
  y    = N - break_after + 1,
  yend = N - break_after + 1,
  x    = break_after,
  xend = N + 0.5
)

# ---- 8) PLOT ----------------------------------------------------------------

# Map short env names to colors, but legend shows long labels with units
env_scale_values <- setNames(
  unname(env_palette_full),
  sub(" \\(.*\\)$", "", names(env_palette_full))
)

env_scale_breaks <- sub(" \\(.*\\)$", "", names(env_palette_full))
env_scale_labels <- names(env_palette_full)

p <- quickcor(node_matrix, type = "upper") +
  
  geom_tile(aes(fill = r)) +
  
  # Region block lines
  geom_segment(
    data = vseg,
    aes(x = x, y = y, xend = xend, yend = yend),
    inherit.aes = FALSE,
    linewidth = 1,
    colour = "black"
  ) +
  geom_segment(
    data = hseg,
    aes(x = x, y = y, xend = xend, yend = yend),
    inherit.aes = FALSE,
    linewidth = 1,
    colour = "black"
  ) +
  
  scale_fill_distiller(
    palette   = "RdBu",
    direction = -1,
    limits    = c(-1, 1),
    breaks    = seq(-1, 1, by = 0.5),
    name      = "Correlation (r)",
    na.value  = "grey90"
  ) +
  
  # Day-only RDA links sized continuously by unique adjusted R² × 100.
  # A heatmap is still saved when no significant positive links qualify.
  (if (nrow(env_links_top) > 0L) {
    anno_link(
      data = env_links_top,
      mapping = aes(colour = .row, size = effect_pct)
    )
  } else {
    NULL
  }) +
  
  scale_colour_manual(
    values = env_scale_values,
    breaks = env_scale_breaks,
    labels = env_scale_labels,
    name   = "Environment"
  ) +
  
  scale_size_continuous(
    limits = c(0, link_size_upper),
    breaks = LINK_LEGEND_BREAKS,
    labels = scales::label_number(suffix = "%"),
    range  = LINK_WIDTH_RANGE,
    name   = "Unique morphological\nvariance explained (%)"
  ) +
  
  guides(
    colour = "none",
    size   = guide_legend(order = 2, override.aes = list(colour = "grey25")),
    fill   = guide_colorbar(order = 3)
  ) +
  
  geom_text(
    data = p_star_df,
    aes(x = x, y = y, label = sig),
    colour = "white",
    fontface = "bold",
    size = 5,
    inherit.aes = FALSE
  ) +
  
  coord_cartesian(clip = "off") +
  
  theme(
    axis.title         = element_blank(),
    axis.text.x.top    = element_blank(),
    axis.text.x.bottom = element_blank(),
    axis.text.y.left   = element_blank(),
    axis.text.y.right  = element_blank(),
    axis.ticks         = element_blank(),
    plot.margin        = margin(120, 180, 30, 180),
    legend.position    = c(0, 0.3)
  ) +
  
  # Region labels, top and right
  annotate(
    "text",
    x = region_blocks_x$xmid,
    y = top_y + 0.9,
    label = region_blocks_x$Region,
    vjust = 0,
    size = 4,
    fontface = "bold"
  ) +
  annotate(
    "text",
    y = region_blocks_y$ymid,
    x = right_x + 0.9,
    label = region_blocks_y$Region,
    angle = 270,
    vjust = 0,
    size = 4,
    fontface = "bold"
  ) +
  
  # Cluster labels, top and right
  annotate(
    "text",
    x = clabs$x,
    y = top_y + 0.2,
    label = clabs$lab,
    vjust = 0,
    size = 4
  ) +
  annotate(
    "text",
    y = ylabs$y,
    x = right_x + 0.2,
    label = ylabs$lab,
    angle = 270,
    vjust = 0,
    size = 4
  )

save_plot_formats(
  file.path(PLOT_DIR, "corplot_envlink_cal_nodeRDA_abundance_weighted.png"),
  p, width = 15, height = 10, dpi = 300
)
