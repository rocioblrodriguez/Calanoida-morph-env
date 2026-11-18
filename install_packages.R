#!/usr/bin/env Rscript
# =============================================================================
# Script:  install_packages.R
# Author:  Rocio Rodriguez
# Purpose: Install dependencies once before running the numbered scripts.
# =============================================================================

# The analysis scripts load their packages without installing or updating them.
packages <- c(
  "readr", "dplyr", "tidyr", "tibble", "purrr", "lubridate", "stringr",
  "sf", "oce", "bestNormalize", "FactoMineR", "vegan", "permute",
  "ggplot2", "ggpubr", "ggtext", "ggimage", "ggpattern", "patchwork", "cowplot",
  "ragg", "svglite", "scales", "Hmisc", "ggcor", "png", "imager", "morphr",
  "rnaturalearth", "rnaturalearthdata", "ggrepel"
)
install.packages(packages, repos = "https://cloud.r-project.org")

# ggcor and morphr may require the repositories used in the original project
# when they are unavailable for the installed R version. See README.md.
