# =============================================================================
# Slope-correlation simulation study
#
# Comparing strategies for estimating the correlation between two longitudinal
# slopes, under complete and (independently) missing data.
# =============================================================================

# -----------------------------------------------------------------------------
# Install and import packages
#
# On a fresh machine you may need system libraries for LMMstar/lava, e.g. on
# Debian/Ubuntu: sudo apt-get install -y libgsl-dev
# Then, run once to install the required R packages:
#   install.packages(c("mvtnorm", "data.table", "boot", "LMMstar", "dplyr", "lava"))
# -----------------------------------------------------------------------------

library(mvtnorm)
library(data.table)
library(boot)
library(LMMstar)
library(dplyr)

# -----------------------------------------------------------------------------
# Scenario 3 -- correlation structure (estimated from ADNI CN subjects)
# -----------------------------------------------------------------------------

n_measurements <- 5

rho_XX <- matrix(
    c(1.00, 0.91, 0.88, 0.86, 0.81,
      0.91, 1.00, 0.89, 0.86, 0.86,
      0.88, 0.89, 1.00, 0.88, 0.89,
      0.86, 0.86, 0.88, 1.00, 0.91,
      0.81, 0.86, 0.89, 0.91, 1.00),
    nrow = 5, byrow = TRUE
)

rho_YY <- matrix(
    c(1.00, 0.91, 0.89, 0.87, 0.86,
      0.91, 1.00, 0.94, 0.92, 0.92,
      0.89, 0.94, 1.00, 0.91, 0.95,
      0.87, 0.92, 0.91, 1.00, 0.94,
      0.86, 0.92, 0.95, 0.94, 1.00),
    nrow = 5, byrow = TRUE
)

rho_YX <- matrix(
    c(-0.14, -0.11, -0.04, -0.04, -0.03,
      -0.14, -0.11, -0.05, -0.04, -0.06,
      -0.14, -0.11, -0.06, -0.07, -0.07,
      -0.12, -0.10, -0.03, -0.06, -0.04,
      -0.09, -0.06, -0.01, -0.03, -0.03),
    nrow = 5, byrow = TRUE
)
rho_XY <- t(rho_YX)

rho <- rbind(cbind(rho_XX, rho_XY),
             cbind(t(rho_XY), rho_YY))

sigma_joint_XY <- diag(10) %*% rho %*% diag(10)

mu <- c(-3.02, -3.31, -3.53, -3.24, -4.16,  # mu_X (PAD)
         41.9,  41.5,  41.2,  40.4,  40.9)  # mu_Y (GM)

# -----------------------------------------------------------------------------
# Contrast matrix -- OLS slope weights for equally-spaced time 0-4
# -----------------------------------------------------------------------------

time   <- 0:(n_measurements - 1)
weight <- (time - mean(time)) / sum((time - mean(time))^2)

contrast <- rbind(c(weight, rep(0, n_measurements)),
                  c(rep(0, n_measurements), weight))

# -----------------------------------------------------------------------------
# True (population) slope correlation
# -----------------------------------------------------------------------------

sigma_bXbY <- contrast %*% sigma_joint_XY %*% t(contrast)
corr_true   <- sigma_bXbY[1, 2] / sqrt(sigma_bXbY[1, 1] * sigma_bXbY[2, 2])
cat("True slope correlation:", round(corr_true, 4), "\n")
# True slope correlation: -0.1695

# -----------------------------------------------------------------------------
# Simulate one dataset and apply independent dropout
# -----------------------------------------------------------------------------

set.seed(42)
n_subjects      <- 10000
dropout_pattern <- c(666, 518, 453, 290, 197)

# Wide format: columns X_1..X_5, Y_1..Y_5
data_wide <- data.table(id = 1:n_subjects,
                        rmvnorm(n_subjects, mean = mu, sigma = sigma_joint_XY))
names(data_wide) <- c("id",
                      paste0("X_", 1:n_measurements),
                      paste0("Y_", 1:n_measurements))

# Long format (complete)
data_long <- melt(data_wide,
                  id.vars       = "id",
                  measure.vars  = patterns("^X_", "^Y_"),
                  value.name    = c("X", "Y"),
                  variable.name = "time")
data_long[, time := as.integer(time) - 1]

# Introduce independent dropout
probs <- c(-diff(dropout_pattern), tail(dropout_pattern, 1)) / dropout_pattern[1]
data_missing <- copy(data_long)
data_missing <- data_missing[, {
    dropout_time <- sample(0:4, size = 1, prob = probs)
    list(time = time,
         X    = ifelse(time > dropout_time, NA, X),
         Y    = ifelse(time > dropout_time, NA, Y))
}, by = id]

cat("\nObservations per time point (with missing data):\n")
print(data_missing[, .(n = sum(!is.na(X))), by = time][order(time)])

# -----------------------------------------------------------------------------
# Strategy 1 -- correlate per-subject OLS slopes
# -----------------------------------------------------------------------------

# -- complete data --
slopes_complete <- data_long %>%
    dplyr::group_by(id) %>%
    dplyr::arrange(time) %>%
    dplyr::summarise(
        slope_x = coef(lm(X ~ time))["time"],
        slope_y = coef(lm(Y ~ time))["time"]
    )
cor_s1_complete <- cor.test(slopes_complete$slope_x, slopes_complete$slope_y)
result_s1_complete <- data.table(
    estimate = cor_s1_complete$estimate,
    lower    = cor_s1_complete$conf.int[1],
    upper    = cor_s1_complete$conf.int[2],
    p.value  = cor_s1_complete$p.value
)
print("Strategy 1, Complete Case")
print(result_s1_complete)

# -- missing data --
slopes_missing <- data_missing %>%
    dplyr::group_by(id) %>%
    dplyr::arrange(time) %>%
    dplyr::summarise(
        slope_x = coef(lm(X ~ time))["time"],
        slope_y = coef(lm(Y ~ time))["time"]
    )
cor_s1_missing <- cor.test(slopes_missing$slope_x, slopes_missing$slope_y)
result_s1_missing <- data.table(
    estimate = cor_s1_missing$estimate,
    lower    = cor_s1_missing$conf.int[1],
    upper    = cor_s1_missing$conf.int[2],
    p.value  = cor_s1_missing$p.value
)
print("Strategy 1, Missing Data")
print(result_s1_missing)

# -----------------------------------------------------------------------------
# Strategy 2a -- sample covariance of all time-point columns + bootstrap CIs
# -----------------------------------------------------------------------------

# Pivot to wide with one column per (variable, time) combination
data_longlong_complete <- melt(data_long,
                               id.vars      = c("id", "time"),
                               measure.vars = c("X", "Y"),
                               variable.name = "variable",
                               value.name    = "value")
data_longlong_complete[, variable := paste(as.character(variable), time, sep = "_")]

data_wide_s2a_complete <- dcast(data_longlong_complete[, -"time"],
                                id ~ variable, value.var = "value")
data_wide_s2a_complete[, id := NULL]

bs_complete <- boot(
    data      = data_wide_s2a_complete,
    stype     = "i",
    statistic = function(d, idx, ctr) {
        sig <- cov(d[idx, ], use = "complete.obs")
        cov2cor(ctr %*% sig %*% t(ctr))[1, 2]
    },
    R        = 1000,
    ctr      = contrast
)
sd_bs_complete  <- sd(bs_complete$t, na.rm = TRUE)
orig_complete   <- bs_complete$t0
result_s2a_complete <- data.table(
    estimate = orig_complete,
    lower    = orig_complete - 1.96 * sd_bs_complete,
    upper    = orig_complete + 1.96 * sd_bs_complete,
    p.value  = 2 * (1 - pnorm(abs(orig_complete / sd_bs_complete)))
)
print("Strategy 2a, Complete Case")
print(result_s2a_complete)

data_longlong_missing <- melt(data_missing,
                              id.vars      = c("id", "time"),
                              measure.vars = c("X", "Y"),
                              variable.name = "variable",
                              value.name    = "value")
data_longlong_missing[, variable := paste(as.character(variable), time, sep = "_")]

data_wide_s2a_missing <- dcast(data_longlong_missing[, -"time"],
                               id ~ variable, value.var = "value")
data_wide_s2a_missing[, id := NULL]

bs_missing <- boot(
    data      = data_wide_s2a_missing,
    stype     = "i",
    statistic = function(d, idx, ctr) {
        sig <- cov(d[idx, ], use = "complete.obs")
        cov2cor(ctr %*% sig %*% t(ctr))[1, 2]
    },
    R        = 1000,
    ctr      = contrast
)
sd_bs_missing  <- sd(bs_missing$t, na.rm = TRUE)
orig_missing   <- bs_missing$t0
result_s2a_missing <- data.table(
    estimate = orig_missing,
    lower    = orig_missing - 1.96 * sd_bs_missing,
    upper    = orig_missing + 1.96 * sd_bs_missing,
    p.value  = 2 * (1 - pnorm(abs(orig_missing / sd_bs_missing)))
)
print("Strategy 2a, Missing Data")
print(result_s2a_missing)

# -----------------------------------------------------------------------------
# Strategy 2b -- covariance via mixed model + delta method CIs
# -----------------------------------------------------------------------------

# -- complete data --
fit_complete <- lmm(value ~ variable,
                    repetition = ~ variable | id,
                    structure  = "UN",
                    data       = data_longlong_complete)
result_s2b_complete <- data.table(
    lava::estimate(fit_complete, function(p) {
        sig       <- sigma(fit_complete, p = p)
        sig_bxby  <- contrast %*% sig %*% t(contrast)
        cov2cor(sig_bxby)[1, 2]
    })
)[, .(estimate, lower, upper, p.value)]

print("Strategy 2b, Complete Case")
print(result_s2b_complete)

# -- missing data --
fit_missing <- lmm(value ~ variable,
                   repetition = ~ variable | id,
                   structure  = "UN",
                   data       = data_longlong_missing)
result_s2b_missing <- data.table(
    lava::estimate(fit_missing, function(p) {
        sig       <- sigma(fit_missing, p = p)
        sig_bxby  <- contrast %*% sig %*% t(contrast)
        cov2cor(sig_bxby)[1, 2]
    })
)[, .(estimate, lower, upper, p.value)]

print("Strategy 2b, Missing Data")
print(result_s2b_missing)

# -----------------------------------------------------------------------------
# Summary
#
# Comparison of the true correlation used for the simulation with estimation
# for the three strategies in presence and absence of missing data.
# -----------------------------------------------------------------------------

all_results <- rbind(result_s1_complete, result_s1_missing,
                     result_s2a_complete, result_s2a_missing,
                     result_s2b_complete, result_s2b_missing)
labels <- c('S1 complete', 'S1 missing', 'S2a complete', 'S2a missing', 'S2b complete', 'S2b missing')
all_results <- cbind(labels, all_results)

print("True correlation:")
print(corr_true)

print("Estimated correlation:")
print(all_results)
