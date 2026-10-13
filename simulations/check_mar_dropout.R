source("R/functions.R")

# Sanity checks for simulate_dropout_mar(). Prints the checks and writes the strategy
# estimates to results/mar_dropout.tsv (plotted by figures/mar_dropout_vis.R).
# Run from the repository root: Rscript simulations/check_mar_dropout.R

n_check <- 200          # datasets for the cheap checks (retention, selection)
n_runs  <- 20           # runs per beta for the strategy comparison (slow: bootstrap + mmrm)
betas   <- c(0, -1, -2) # beta < 0: subjects losing more GM (Y) drop out more

n_subjects     <- 1000
n_measurements <- 5
dropout_pattern <- c(666, 518, 453, 290, 197) # taken from real data, adjusted to be monotonic

time   <- 0:(n_measurements - 1)
weight <- (time - mean(time)) / sum((time - mean(time))^2)
contrast <- rbind(c(weight, rep(0, n_measurements)), c(rep(0, n_measurements), weight))

# Scenario 3 setup (copied from simulations/simulated_data.R)
sigma_marginalXY <- diag(rep(1, 2 * n_measurements))

rho_YY <- matrix(
    c(
        1.00, 0.91, 0.89, 0.87, 0.86,
        0.91, 1.00, 0.94, 0.92, 0.92,
        0.89, 0.94, 1.00, 0.91, 0.95,
        0.87, 0.92, 0.91, 1.00, 0.94,
        0.86, 0.92, 0.95, 0.94, 1.00
    ),
    nrow = 5, ncol = 5, byrow = TRUE
)

rho_XX <- matrix(
    c(
        1.00, 0.91, 0.88, 0.86, 0.81,
        0.91, 1.00, 0.89, 0.86, 0.86,
        0.88, 0.89, 1.00, 0.88, 0.89,
        0.86, 0.86, 0.88, 1.00, 0.91,
        0.81, 0.86, 0.89, 0.91, 1.00
    ),
    nrow = 5, ncol = 5, byrow = TRUE
)

rho_YX <- matrix(
    c(
        -0.14, -0.11, -0.04, -0.04, -0.03,
        -0.14, -0.11, -0.05, -0.04, -0.06,
        -0.14, -0.11, -0.06, -0.07, -0.07,
        -0.12, -0.10, -0.03, -0.06, -0.04,
        -0.09, -0.06, -0.01, -0.03, -0.03
    ),
    nrow = 5, ncol = 5, byrow = TRUE
)

rho_XY <- t(rho_YX)
rho <- rbind(cbind(rho_XX, rho_XY), cbind(t(rho_XY), rho_YY))

mu <- c(-3.02, -3.31, -3.53, -3.24, -4.16,   # mu_X (PAD)
         41.9,  41.5,  41.2,  40.4,  40.9)    # mu_Y (GM)

sigma_joint_XY       <- sigma_marginalXY %*% rho %*% sigma_marginalXY
sigma_bXbY_realistic <- contrast %*% sigma_joint_XY %*% t(contrast)
corr_true            <- cov2cor(sigma_bXbY_realistic)[1, 2]

simulate_one <- function(seed, beta) {
    set.seed(seed)
    data_long <- make_data_long(simulate_data(n_subjects, mu, sigma_joint_XY, n_measurements))
    list(full = data_long, missing = simulate_dropout_mar(data_long, dropout_pattern, beta))
}

# Check 1: realized retention per visit matches dropout_pattern for every beta.
message("\n== Check 1: retention per visit (mean number observed over ", n_check, " datasets) ==")
retention <- rbindlist(lapply(betas, function(beta) {
    counts <- sapply(1:n_check, function(i) {
        d <- simulate_one(i, beta)$missing
        d[order(time), sum(!is.na(Y)), by = time]$V1
    })
    data.table(beta = beta, time = time, observed = rowMeans(counts))
}))
retention <- dcast(retention, time ~ paste0("beta_", beta), value.var = "observed")
retention[, expected := dropout_pattern / dropout_pattern[1] * n_subjects]
print(retention, digits = 4)

# Check 2: selection. With beta < 0, completers should have lost less Y than everyone.
message("\n== Check 2: mean change Y_4 - Y_0, all subjects vs completers ==")
selection <- rbindlist(lapply(betas, function(beta) {
    rbindlist(lapply(1:n_check, function(i) {
        d <- simulate_one(i, beta)
        change <- d$full[order(id, time), Y[.N] - Y[1], by = id]
        completers <- d$missing[time == max(time) & !is.na(Y), id]
        data.table(all = mean(change$V1), completers = mean(change[id %in% completers, V1]))
    }))[, .(beta = beta, all = mean(all), completers = mean(completers))]
}))
print(selection, digits = 4)

# Checks 3 and 4: strategy bias per beta. beta = 0 should match the existing MCAR
# dropout; as |beta| grows, Strategy 2b (mmrm, valid under MAR) should stay near the
# truth while Strategies 1 and 2a may drift. If 2b drifts, suspect an MNAR leak.
message("\n== Checks 3 + 4: strategy estimates over ", n_runs, " runs (true correlation = ", round(corr_true, 4), ") ==")
runs <- c(
    list(run_simulation(n_runs, n_subjects, n_measurements, mu, sigma_joint_XY, contrast,
                        dropout_pattern = dropout_pattern, mechanism = "dropout",
                        label = "MCAR dropout (reference)")[, mechanism := "MCAR"]),
    lapply(betas, function(beta) {
        run_simulation(n_runs, n_subjects, n_measurements, mu, sigma_joint_XY, contrast,
                       dropout_pattern = dropout_pattern, mechanism = "dropout_mar", beta_mar = beta,
                       label = paste0("MAR dropout, beta = ", beta))[, mechanism := paste0("MAR beta=", beta)]
    })
)
estimates <- rbindlist(runs)
estimates[, value := as.numeric(value)]
# set estimates outside [-1, 1] to NA (numerical instability in small bootstrap resamples)
estimates[value < -1 | value > 1, value := NA]
summary_table <- estimates[, .(mean = mean(value, na.rm = TRUE),
                               bias = mean(value, na.rm = TRUE) - corr_true,
                               sd   = sd(value, na.rm = TRUE)),
                           by = .(mechanism, strategy = variable)]
print(dcast(summary_table, strategy ~ mechanism, value.var = "bias"), digits = 3)
message("(values above are bias; full table below)")
print(summary_table, digits = 3)

write.table(estimates[, .(experiment = mechanism, strategy = variable, value, true = corr_true)],
            "results/mar_dropout.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
print(summary(warnings()))
