# =============================================================================
# run_all.R — reproduce all simulation results and figures
#
# Run from the project root:
#   Rscript run_all.R
# or interactively in RStudio with the project open (working dir = project root).
#
# Output: all PDFs and TSVs written to results/
#
# Runtime: scales with n_runs in simulations/simulated_data.R (6 scenarios × 2 adjustments × 3 strategies).
# For a quick smoke-test, open simulations/simulated_data.R and temporarily
# set n_runs <- 5 before sourcing.
# =============================================================================

if (!dir.exists("results")) dir.create("results")

# ── Simulations ───────────────────────────────────────────────────────────────
# Each script sources R/functions.R and writes to results/.
# simulated_data.R must run before simulation_vis.R (produces the TSV).

message("=== Running simulation: main scenarios ===")
source("simulations/simulated_data.R")

message("=== Performance tables ===")
source("simulations/performance_table.R")

# ── Figures ───────────────────────────────────────────────────────────────────

message("=== Figure: simulation results ===")
source("figures/simulation_vis.R")

message("=== Figure: covariance heatmaps ===")
source("figures/heatmap_sigma.R")

message("=== Done. All outputs are in results/ ===")
