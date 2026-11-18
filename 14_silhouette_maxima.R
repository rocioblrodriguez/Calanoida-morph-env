#!/usr/bin/env Rscript
# =============================================================================
# Script:  14_silhouette_maxima.R
# Author:  Rocio Rodriguez
# Purpose: Summarize global and secondary local silhouette maxima.
# Input:   Morphology diagnostic CSV at DIAGNOSTICS_CSV (00_setup.R).
# Outputs: Maxima/frequency tables and silhouette-profile figures.
# =============================================================================

# 1. Settings and packages -----------------------------------------------------
source("R_Scripts/00_setup.R", local = TRUE)
PLOT_DIR <- PLOT_DIRS[["Silhouette"]]
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# 2. Read archived diagnostic results -----------------------------------------
# Each row describes one candidate k within one variable-sampling iteration.
# This script summarizes those results; it does not regenerate the upstream
# variable-sampling/k-means diagnostics or choose the final 20 descriptors.
dat <- read_csv(DIAGNOSTICS_CSV, show_col_types = FALSE) |>
  select(iteration, k, silhouette_avg) |>
  arrange(iteration, k)
n_iterations <- n_distinct(dat$iteration)

# 3. Identify the global maximum in each iteration -----------------------------
# Ties select the smallest k because the input is ordered by k.
first_max <- dat |>
  group_by(iteration) |>
  slice_max(silhouette_avg, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(maximum_type = "First maximum")

# 4. Identify the first internal local maximum after the global maximum --------
# A true internal peak exceeds both adjacent k values. Endpoints have only one
# neighbor and cannot be secondary peaks. Later internal peaks are not counted.
second_max <- dat |>
  left_join(select(first_max, iteration, global_k = k), by = "iteration") |>
  group_by(iteration) |>
  mutate(is_local_max = silhouette_avg > lag(silhouette_avg) &
                         silhouette_avg > lead(silhouette_avg)) |>
  filter(k > global_k, is_local_max) |>
  slice_min(k, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(iteration, k, silhouette_avg) |>
  mutate(maximum_type = "Second maximum")

# 5. Count maxima and report their occurrence ---------------------------------
counts <- bind_rows(first_max, second_max) |>
  count(k, maximum_type, name = "count") |>
  complete(k = 2:10, maximum_type = c("First maximum", "Second maximum"), fill = list(count = 0L)) |>
  mutate(maximum_type = factor(maximum_type, levels = c("First maximum", "Second maximum")))
write_csv(first_max, file.path(OUT_DIR, "silhouette_FIRST_global_maximum.csv"))
write_csv(second_max, file.path(OUT_DIR, "silhouette_SECOND_maximum_by_iteration.csv"))
write_csv(counts, file.path(OUT_DIR, "silhouette_first_second_maximum_counts.csv"))
write_csv(filter(counts, maximum_type == "Second maximum"),
          file.path(OUT_DIR, "silhouette_SECOND_maximum_counts.csv"))
message("Iterations: ", n_iterations, "; iterations with a secondary maximum: ", nrow(second_max))

# 6. Original frequency panel (display k = 2–9, as supplied) ------------------
plot_dat <- filter(counts, k %in% 2:9)
p <- ggplot(
  plot_dat,
  aes(
    x = factor(k),
    y = count,
    fill = maximum_type
  )
) +
  
  geom_col(
    width = 0.70,
    colour = "white",
    linewidth = 0.3
  ) +
  
  scale_fill_manual(
    
    values = c(
      
      # k=2 global maximum
      "First maximum" =
        "#D7191C",
      
      # first internal local maximum after k=2
      "Second maximum" =
        "#2C7BB6"
    ),
    
    labels = c(
      "First maximum",
      "Second maximum"
    ),
    
    name = NULL
  ) +
  
  scale_x_discrete(
    drop = FALSE
  ) +
  
  scale_y_continuous(
    expand = expansion(
      mult = c(0, 0.08)
    )
  ) +
  
  labs(
    x = "Number of clusters (k)",
    y = "Number of iterations"
  ) +
  
  theme_minimal(
    base_size = 13
  ) +
  
  theme(
    
    legend.position = "top",
    
    panel.grid.major.x =
      element_blank(),
    
    axis.ticks.x =
      element_blank()
  )


save_plot_formats(file.path(PLOT_DIR, "silhouette_first_second_maximum.png"),
                  p, width = 8, height = 5, dpi = 300)

# 7. Original profile panel ---------------------------------------------------
dat_k6_second <- semi_join(dat, filter(second_max, k == 6) |> select(iteration), by = "iteration")
iterations_k6_second <- unique(dat_k6_second$iteration)
mean_k6_second <- dat_k6_second |>
  group_by(k) |>
  summarise(mean_silhouette = mean(silhouette_avg), .groups = "drop")
p6 <- ggplot() +
  geom_line(
    data = dat_k6_second,
    aes(
      x = k,
      y = silhouette_avg,
      group = iteration
    ),
    colour = "#2C7BB6",
    alpha = 0.15,
    linewidth = 0.4
  ) +
  geom_line(
    data = mean_k6_second,
    aes(
      x = k,
      y = mean_silhouette
    ),
    colour = "black",
    linewidth = 1.4
  ) +
  geom_point(
    data = mean_k6_second,
    aes(
      x = k,
      y = mean_silhouette
    ),
    colour = "black",
    size = 2.8
  ) +
  geom_vline(
    xintercept = 6,
    linetype = "dashed",
    linewidth = 0.7,
    colour = "black"
  ) +
  scale_x_continuous(
    breaks = 2:10
  ) +
  labs(
    subtitle = paste0(
      "k = 6 secondary local maximum; n = ",
      length(iterations_k6_second),
      " iterations"
    ),
    x = "Number of clusters (k)",
    y = "Mean silhouette width"
  ) +
  theme_classic(
    base_size = 13
  )

save_plot_formats(file.path(PLOT_DIR, "silhouette_profiles_k6_SECOND_maximum.png"),
                  p6, width = 8, height = 5, dpi = 300)
write_csv(mean_k6_second, file.path(OUT_DIR, "silhouette_profiles_k6_SECOND_maximum_mean.csv"))
