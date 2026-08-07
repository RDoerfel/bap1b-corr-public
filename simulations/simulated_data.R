source("R/functions.R")
set.seed(42)

n_subjects    <- 1000
n_measurements <- 5
n_runs        <- 100
dropout_pattern <- c(666, 518, 453, 290, 197) # taken from real data, adjusted to be monotonic

time   <- 0:(n_measurements - 1)
weight <- (time - mean(time)) / sum((time - mean(time))^2)
contrast <- rbind(c(weight, rep(0, n_measurements)), c(rep(0, n_measurements), weight))
contrast_adjusted <- rbind(
    c(weight, rep(0, n_measurements)),
    c(rep(0, n_measurements), weight),
    c(1, rep(0, n_measurements - 1), rep(0, n_measurements)),
    c(rep(0, n_measurements), c(1, rep(0, n_measurements - 1)))
)

mu_X <- 3 + 0:(n_measurements - 1)
mu_Y <- 5 - 0:(n_measurements - 1)
mu   <- c(mu_X, mu_Y)

sigma_marginalXY <- diag(rep(1, 2 * n_measurements))

# Scenario 1: constant diagonal cross-correlation
rho_X          <- diag(1, n_measurements)
rho_Y          <- diag(1, n_measurements)
rho_XY_diagonal <- diag(0.6, n_measurements)

rho_diag          <- rbind(cbind(rho_X, t(rho_XY_diagonal)), cbind(rho_XY_diagonal, rho_Y))
sigma_XY_diag_const <- sigma_marginalXY %*% rho_diag %*% sigma_marginalXY

cov_true_bxby_diag_const <- contrast %*% sigma_XY_diag_const %*% t(contrast)
corr_true_diag_const <- cov_true_bxby_diag_const[1, 2] / sqrt(cov_true_bxby_diag_const[1, 1] * cov_true_bxby_diag_const[2, 2])

cov_true_bxby_diag_const_adjusted <- contrast_adjusted %*% sigma_XY_diag_const %*% t(contrast_adjusted)
cov_true_bxby_diag_const_adjusted_final <- cov_true_bxby_diag_const_adjusted[1:2, 1:2] - # sigma_slope
                                     cov_true_bxby_diag_const_adjusted[1:2, 3:4] %*%  # sigma_cross
                                     solve(cov_true_bxby_diag_const_adjusted[3:4, 3:4]) %*% # sigma_baseline^(-1)
                                     cov_true_bxby_diag_const_adjusted[3:4, 1:2] # sigma_cross^T
corr_true_diag_const_adjusted <- cov_true_bxby_diag_const_adjusted_final[1, 2] / sqrt(cov_true_bxby_diag_const_adjusted_final[1, 1] * cov_true_bxby_diag_const_adjusted_final[2, 2])

results_diag_const <- run_simulation(n_runs, n_subjects, n_measurements, mu, sigma_XY_diag_const, contrast, dropout_pattern = dropout_pattern, mechanism = "dropout", label = "Scenario 1 (unadjusted)")
results_diag_const_adjusted <- run_simulation(n_runs, n_subjects, n_measurements, mu, sigma_XY_diag_const, contrast_adjusted, adjusted = TRUE, dropout_pattern = dropout_pattern, mechanism = "dropout", label = "Scenario 1 (adjusted)")

# Scenario 2: non-constant diagonal cross-correlation
rho_XY_diagonal <- diag(c(0.8, 0.7, 0.6, 0.5, 0.4))

rho_diag     <- rbind(cbind(rho_X, t(rho_XY_diagonal)), cbind(rho_XY_diagonal, rho_Y))
sigma_XY_diag <- sigma_marginalXY %*% rho_diag %*% sigma_marginalXY

cov_true_bxby_diag <- contrast %*% sigma_XY_diag %*% t(contrast)
corr_true_diag <- cov_true_bxby_diag[1, 2] / sqrt(cov_true_bxby_diag[1, 1] * cov_true_bxby_diag[2, 2])

cov_true_bxby_diag_adjusted <- contrast_adjusted %*% sigma_XY_diag %*% t(contrast_adjusted)
cov_true_bxby_diag_adjusted_final <- cov_true_bxby_diag_adjusted[1:2, 1:2] - # sigma_slope
                                     cov_true_bxby_diag_adjusted[1:2, 3:4] %*%  # sigma_cross
                                     solve(cov_true_bxby_diag_adjusted[3:4, 3:4]) %*% # sigma_baseline^(-1)
                                     cov_true_bxby_diag_adjusted[3:4, 1:2] # sigma_cross^T
corr_true_diag_adjusted <- cov_true_bxby_diag_adjusted_final[1, 2] / sqrt(cov_true_bxby_diag_adjusted_final[1, 1] * cov_true_bxby_diag_adjusted_final[2, 2])

results_diag          <- run_simulation(n_runs, n_subjects, n_measurements, mu, sigma_XY_diag, contrast, dropout_pattern = dropout_pattern, mechanism = "dropout", label = "Scenario 2 (unadjusted)")
results_diag_adjusted <- run_simulation(n_runs, n_subjects, n_measurements, mu, sigma_XY_diag, contrast_adjusted, adjusted = TRUE, dropout_pattern = dropout_pattern, mechanism = "dropout", label = "Scenario 2 (adjusted)")

# Scenario 3: realistic covariance structure, estimated from ADNI (CN subjects)
rho_YY <- matrix(
    c(
        1.00, 0.91, 0.89, 0.87, 0.86,
        0.91, 1.00, 0.94, 0.92, 0.92,
        0.89, 0.94, 1.00, 0.91, 0.95,
        0.87, 0.92, 0.91, 1.00, 0.94,
        0.86, 0.92, 0.95, 0.94, 1.00
    ),
    nrow = 5, ncol = 5, byrow = TRUE,
    dimnames = list(paste0("Y_", 0:4), paste0("Y_", 0:4))
)

rho_XX <- matrix(
    c(
        1.00, 0.91, 0.88, 0.86, 0.81,
        0.91, 1.00, 0.89, 0.86, 0.86,
        0.88, 0.89, 1.00, 0.88, 0.89,
        0.86, 0.86, 0.88, 1.00, 0.91,
        0.81, 0.86, 0.89, 0.91, 1.00
    ),
    nrow = 5, ncol = 5, byrow = TRUE,
    dimnames = list(paste0("X_", 0:4), paste0("X_", 0:4))
)

# Cross-correlation block: Y (gm) rows, X (pad) columns
rho_YX <- matrix(
    c(
        -0.14, -0.11, -0.04, -0.04, -0.03,
        -0.14, -0.11, -0.05, -0.04, -0.06,
        -0.14, -0.11, -0.06, -0.07, -0.07,
        -0.12, -0.10, -0.03, -0.06, -0.04,
        -0.09, -0.06, -0.01, -0.03, -0.03
    ),
    nrow = 5, ncol = 5, byrow = TRUE,
    dimnames = list(paste0("Y_", 0:4), paste0("X_", 0:4))
)

rho_XY <- t(rho_YX)
rho <- rbind(cbind(rho_XX, rho_XY), cbind(t(rho_XY), rho_YY))

mu <- c(-3.02, -3.31, -3.53, -3.24, -4.16,   # mu_X (PAD)
         41.9,  41.5,  41.2,  40.4,  40.9)    # mu_Y (GM)

sigma_joint_XY    <- sigma_marginalXY %*% rho %*% sigma_marginalXY
sigma_bXbY_realistic <- contrast %*% sigma_joint_XY %*% t(contrast)
corr_true_realistic  <- sigma_bXbY_realistic[1, 2] / sqrt(sigma_bXbY_realistic[1, 1] * sigma_bXbY_realistic[2, 2])

sigma_bXbY_realistic_adjusted <- contrast_adjusted %*% sigma_joint_XY %*% t(contrast_adjusted)
sigma_bXbY_realistic_adjusted_final <- sigma_bXbY_realistic_adjusted[1:2, 1:2] - # sigma_slope
                                     sigma_bXbY_realistic_adjusted[1:2, 3:4] %*%  # sigma_cross
                                     solve(sigma_bXbY_realistic_adjusted[3:4, 3:4]) %*% # sigma_baseline^(-1)
                                     sigma_bXbY_realistic_adjusted[3:4, 1:2] # sigma_cross^T
corr_true_realistic_adjusted <- sigma_bXbY_realistic_adjusted_final[1, 2] / sqrt(sigma_bXbY_realistic_adjusted_final[1, 1] * sigma_bXbY_realistic_adjusted_final[2, 2])

results_realistic          <- run_simulation(n_runs, n_subjects, n_measurements, mu, sigma_joint_XY, contrast, dropout_pattern = dropout_pattern, mechanism = "dropout", label = "Scenario 3 (unadjusted)")
results_realistic_adjusted <- run_simulation(n_runs, n_subjects, n_measurements, mu, sigma_joint_XY, contrast_adjusted, adjusted = TRUE, dropout_pattern = dropout_pattern, mechanism = "dropout", label = "Scenario 3 (adjusted)")

results <- rbind(
    data.table(experiment = "Scenario 1", adjusted = FALSE, strategy = results_diag_const$variable,          value = results_diag_const$value,          true = corr_true_diag_const),
    data.table(experiment = "Scenario 2", adjusted = FALSE, strategy = results_diag$variable,                value = results_diag$value,                true = corr_true_diag),
    data.table(experiment = "Scenario 3", adjusted = FALSE, strategy = results_realistic$variable,           value = results_realistic$value,           true = corr_true_realistic),
    data.table(experiment = "Scenario 1", adjusted = TRUE,  strategy = results_diag_const_adjusted$variable, value = results_diag_const_adjusted$value,  true = corr_true_diag_const_adjusted),
    data.table(experiment = "Scenario 2", adjusted = TRUE,  strategy = results_diag_adjusted$variable,       value = results_diag_adjusted$value,        true = corr_true_diag_adjusted),
    data.table(experiment = "Scenario 3", adjusted = TRUE,  strategy = results_realistic_adjusted$variable,  value = results_realistic_adjusted$value,   true = corr_true_realistic_adjusted)
)

results$value <- as.numeric(results$value)
# set estimates outside [-1, 1] to NA (numerical instability in small bootstrap resamples)
results$value[results$value < -1 | results$value > 1] <- NA

write.table(results, "results/correlation_slope_vs_slope.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
