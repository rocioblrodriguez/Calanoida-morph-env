# CHP1 revised workflow — Subarctic ST excluded

Run from this project directory with:

```r
source("R_Scripts/00_run_all.R")
```

The six-cluster model is fitted and labelled using all validated NH1208
Calanoida. After labels are assigned, observations with `Region ==
"Subarctic"` and `Cluster == "ST"` are removed before the clustered dataset
and downstream results are saved. Consequently, the PCA, size/transparency,
behavior, RDA, centroid/environment-link, butterfly, and random-image outputs
all exclude Subarctic ST. ST remains in Transition and Subtropical.

The morphology-only correlation diagnostic in script 04 reads the filtered
clustered dataset. Sampling maps and environmental sections are not
cluster-based and therefore retain their full sampling data.

Input paths and shared settings are defined in `R_Scripts/00_setup.R`.
Figures are exported as PNG, PDF, SVG, and EPS.
