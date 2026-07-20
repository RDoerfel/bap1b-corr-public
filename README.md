# Correlation of Rates of Change in Longitudinal Studies with Dropout

Code accompanying the paper:

> **Estimating the correlation between rates of change in longitudinal studies in the presence of dropout**  
> *[Authors] — [Journal, Year] — DOI: [TODO]*

## Overview

The common *slopes-first* estimator (fit per-subject OLS regressions, then correlate the resulting slopes) is biased under independent dropout whenever the cross-covariance between two longitudinal signals is not constant over time. This repository provides:

- The accompanying code for the simulation results of the paper.
- A tutorial for users who want to apply the estimators to their own data.

## Repository structure

```
bap1b-corr-public/
├── R/                        # core statistical functions (sourced by all scripts)
│   └── functions.R
├── simulations/              # Monte Carlo simulation scripts
│   ├── simulated_data.R      # main bias study: Scenarios 1–3
│   └── simulated_adjusted.R  # baseline-adjusted correlation (single dataset)
├── figures/                  # figure-generation scripts (run after simulations)
│   ├── simulation_vis.R      # Fig: simulation results (reads results/…tsv)
│   └── heatmap_sigma.R       # Fig: covariance structure heatmaps
├── real_data/
│   └── real_data.R           # real ADNI analysis (illustrative; data not shared)
├── tutorial/
│   ├── tutorial_script.R     # self-contained R tutorial
│   └── tutorial_colab.ipynb  # Google Colab notebook version
├── results/                  # generated outputs (gitignored; reproduced by run_all.R)
├── run_all.R                 # single entry point to reproduce everything
├── INSTALL.R                 # installs required R packages
└── LICENSE
```

## Tutorial

`tutorial/tutorial_script.R` walks through Scenario 3 (realistic ADNI-derived covariance structure) step by step, comparing all three estimators on both complete and dropout data.

An identical notebook version can be run directly in the browser:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/RDoerfel/bap1b-corr-public/blob/main/tutorial/tutorial_colab.ipynb)

## Reproducing the results

### 1 — System requirements

R ≥ 4.2 is required. On Debian/Ubuntu, install the system libraries needed by `LMMstar` and `mmrm`:

```bash
sudo apt-get install -y libgsl-dev liblapack-dev libblas-dev
```

### 2 — Install R packages

```r
Rscript INSTALL.R
```

### 3 — Run all simulations and generate all figures

```r
Rscript run_all.R
```

All outputs (TSV tables and PDF figures) are written to `results/`.

**Runtime:** approximately 30–60 minutes (100 Monte Carlo runs × 3 scenarios × 2 variants). For a quick smoke-test, temporarily set `n_runs <- 5` inside `simulations/simulated_data.R`.

### Running scripts individually

Each script assumes the **project root** as the working directory. In RStudio, open the `.Rproj` file (or set the working directory manually). Scripts can then be sourced individually in the order listed in `run_all.R`.

## Real data analysis

`real_data/real_data.R` documents the full analysis pipeline applied to the ADNI dataset (PAD vs. grey matter volume in CN subjects). The underlying data are not publicly available and cannot be shared; the script is provided for transparency so that the reported results can be understood and verified against the methodology in the paper.

## License

MIT — see [LICENSE](LICENSE).
