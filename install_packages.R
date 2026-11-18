#!/usr/bin/env Rscript
# =============================================================================
# Script:  install_packages.R
# Author:  Rocio Rodriguez
# Purpose: Install dependencies once before running the numbered scripts.
# =============================================================================

# -----------------------------------------------------------------------------
# Workflow notes:
# - Installs the R packages required by the Chapter 1 workflow before the numbered scripts are run.
# - Package installation is intentionally separated from analysis execution so routine reruns do not update the software environment.
# - Project-specific packages that may not be available for every R version are noted at the end of this file.
# -----------------------------------------------------------------------------


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
