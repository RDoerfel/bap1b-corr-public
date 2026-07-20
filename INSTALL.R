# =============================================================================
# INSTALL.R — install required R packages
#
# Run once before running run_all.R or the tutorial:
#   Rscript INSTALL.R
#
# System dependencies (Debian/Ubuntu):
#   sudo apt-get install -y libgsl-dev liblapack-dev libblas-dev
# =============================================================================

pkgs <- c(
  # simulation study
  "mvtnorm",
  "data.table",
  "boot",
  "LMMstar",    # also installs lava
  "lava",
  "mmrm",
  "dplyr",
  "ggplot2",
  "patchwork",
  "gridExtra",
  "reshape2",
  # real data analysis (requires access to data/data.xlsx)
  "readxl",
  "prodlim",
  "riskRegression",
  "lme4",
  "lmerTest"
)

install.packages(pkgs, repos = "https://cloud.r-project.org")
