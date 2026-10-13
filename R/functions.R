library(mvtnorm)
library(data.table)
library(LMMstar)
library(boot)
library(mmrm)
library(dplyr)

# Simulate data for n_subjects subjects, with a multivariate normal distribution.
simulate_data <- function(n_subjects, mu, sigma, n_measurements, baseline_case=FALSE) {
    data_wide <- data.table(id = 1:n_subjects, rmvnorm(n_subjects, mean = mu, sigma = sigma))
    if (baseline_case) {
        names(data_wide) <- c("id", paste0("X_", 1), paste0("Y_", 1:(n_measurements)))
    } else {
        names(data_wide) <- c("id", paste0("X_", 1:n_measurements), paste0("Y_", 1:n_measurements))
    }
    return(data_wide)
}

# Convert data_wide to data_long format, so that X and Y are in separate columns.
make_data_long <- function(data_wide, baseline_case=FALSE) {
    if (baseline_case) {
        data_long <- melt(data_wide,
                          id.vars = c("id", "X_1"),
                          measure.vars = patterns("^Y_"),
                          value.name = "Y",
                          variable.name = "time")
        setnames(data_long, "X_1", "X")
        data_long[, time := as.integer(time) - 1]
    } else {
        data_long <- melt(data_wide,
                          id.vars = "id",
                          measure.vars = patterns("^X_", "^Y_"),
                          value.name = c("X", "Y"),
                          variable.name = "time")
        data_long[, time := as.integer(time) - 1]
    }
    return(data_long)
}

# Convert data_long to data_longlong format, so that X and Y are in the same column.
make_data_longlong <- function(data_long, baseline_case=FALSE) {
    data_longlong <- melt(data_long,
                          id.vars = c("id", "time"),
                          measure.vars = c("X", "Y"),
                          variable.name = "variable",
                          value.name = "value")
    data_longlong[, variable := as.character(variable)]
    if (baseline_case) {
        data_longlong <- data_longlong[variable != "X" | time == 0]
        data_longlong[, variable := ifelse(variable == "X", "X_0", paste0("Y_", time))]
    } else {
        data_longlong[, variable := paste(variable, time, sep = "_")]
    }
    return(data_longlong)
}

# Simulate dropout based on ADNI data for Healthy subjects.
simulate_dropout <- function(data_long, dropout_pattern) {
    data_long_missing <- copy(data_long)

    # Probability of each timepoint being a subject's last observation
    probs <- c(-diff(dropout_pattern), tail(dropout_pattern, 1)) / dropout_pattern[1]

    data_long_missing[, {
        dropout_time <- sample(0:4, size = 1, prob = probs)
        X_new <- ifelse(time > dropout_time, NA, X)
        Y_new <- ifelse(time > dropout_time, NA, Y)

        list(time = time, X = X_new, Y = Y_new)
    }, by = id]
}

# Simulate monotone dropout that is missing at random (MAR). After each visit t,
# a subject still in the study drops out with probability plogis(a_t + beta * z_t),
# where z_t is the standardized change from baseline Y_t - Y_0 (observed data only,
# so the mechanism is MAR). The intercepts a_t are calibrated per dataset so the
# expected retention matches dropout_pattern; beta = 0 reproduces MCAR dropout.
simulate_dropout_mar <- function(data_long, dropout_pattern, beta) {
    n_visits <- length(dropout_pattern)

    # Conditional probability of dropping out after visit t, given still in at t
    hazard <- 1 - dropout_pattern[-1] / dropout_pattern[-n_visits]

    # Wide matrix of Y (one row per subject, one column per visit) for fast per-visit access
    Y_wide <- dcast(data_long, id ~ time, value.var = "Y")
    ids <- Y_wide$id
    Y_mat <- as.matrix(Y_wide[, -1])

    # Everyone starts as a completer (last visit) and at risk of dropping out
    subject_dropout_time <- rep(n_visits - 1, length(ids))
    at_risk <- rep(TRUE, length(ids))

    # Decide dropout visit by visit, using only data observed up to visit t (keeps it MAR)
    for (t in 0:(n_visits - 2)) {
        if (t == 0) {
            z <- rep(0, length(ids)) # no change at baseline, so the first transition is MCAR
        } else {
            # standardize with the complete data so beta has the same meaning at every visit
            change <- Y_mat[, t + 1] - Y_mat[, 1]
            z <- (change - mean(change)) / sd(change)
        }
        # Nothing to do if nobody should drop at this visit or nobody is left
        if (hazard[t + 1] <= 0 || !any(at_risk)) next

        # Pick the intercept so the average dropout probability equals the ADNI hazard
        z_risk <- z[at_risk]
        intercept <- uniroot(function(a) mean(plogis(a + beta * z_risk)) - hazard[t + 1],
                             interval = c(-50, 50))$root
        # Each at-risk subject drops with a probability that depends on their own change in Y
        drop <- rbinom(length(z_risk), size = 1, prob = plogis(intercept + beta * z_risk)) == 1

        # Record the last observed visit for those who drop and remove them from the risk set
        dropped <- which(at_risk)[drop]
        subject_dropout_time[dropped] <- t
        at_risk[dropped] <- FALSE
    }

    # Set X and Y to NA after each subject's last observed visit (monotone dropout)
    data_long %>%
        dplyr::mutate(
            dropout_time = subject_dropout_time[match(id, ids)],
            X = ifelse(time > dropout_time, NA, X),
            Y = ifelse(time > dropout_time, NA, Y)
        ) %>%
        dplyr::select(-dropout_time)
}

# Simulate intermittent missingness based on the same ADNI-derived pattern,
# without the monotone dropout constraint: each timepoint is independently
# missing with the marginal probability implied by the pattern, so a subject
# can miss a visit and still be observed at a later one.
simulate_missingness <- function(data_long, dropout_pattern) {
    data_long_missing <- copy(data_long)

    # Marginal probability of being missing at each timepoint
    miss_prob <- 1 - dropout_pattern / dropout_pattern[1]

    data_long_missing[, {
        missing <- rbinom(length(time), size = 1, prob = miss_prob[time + 1]) == 1
        X_new <- ifelse(missing, NA, X)
        Y_new <- ifelse(missing, NA, Y)

        list(time = time, X = X_new, Y = Y_new)
    }, by = id]
}

# Strategy 1: correlate per-subject OLS slopes.
strategy_1 <- function(data_long, baseline_case=FALSE) {
    if (baseline_case) {
        slopes <- data_long %>%
            dplyr::group_by(id) %>%
            dplyr::arrange(time) %>%
            dplyr::summarise(
                baseline_x = dplyr::first(X),
                slope_y = coef(lm(Y ~ time))["time"]
            )
        cor <- cor.test(slopes$baseline_x, slopes$slope_y)
    } else {
        slopes <- data_long %>%
            dplyr::group_by(id) %>%
            dplyr::arrange(time) %>%
            dplyr::summarise(
                slope_x = coef(lm(X ~ time))["time"],
                slope_y = coef(lm(Y ~ time))["time"]
            )
        cor <- cor.test(slopes$slope_x, slopes$slope_y)
    }
    table_cor <- data.table(estimate = cor$estimate, lower = cor$conf.int[1], upper = cor$conf.int[2], p.value = cor$p.value)

    return(table_cor)
}

strategy_1_adjusted <- function(data_long) {
    slopes <- data_long %>%
        dplyr::group_by(id) %>%
        dplyr::arrange(time) %>%
        dplyr::summarise(
            slope_x = coef(lm(X ~ time))["time"],
            slope_y = coef(lm(Y ~ time))["time"],
            baseline_x = X[1],
            baseline_y = Y[1],
        )

    lm_adjusted <- lm(formula = slope_y ~ slope_x + baseline_x + baseline_y, data = slopes)
    estimate <- coef(lm_adjusted)["slope_x"]
    lower <- confint(lm_adjusted)["slope_x", 1]
    upper <- confint(lm_adjusted)["slope_x", 2]
    p.value <- summary(lm_adjusted)$coefficients["slope_x", 4]
    table_cor <- data.table(estimate = estimate, lower = lower, upper = upper, p.value = p.value)

    return(table_cor)
}

# Using LMM to estimate the joint variance covariance matrix.
# This is slower than mmrm, but provides confidence intervals for the correlation.
strategy_2b <- function(data_longlong, contrast) {
    unstructured_lmm <- lmm(value ~ variable, repetition = ~ variable | id, structure = "UN", data = data_longlong)
    cor <- lava::estimate(unstructured_lmm, function(p){
        sigma_xy <- sigma(unstructured_lmm, p = p)
        sigma_bxby <- contrast %*% sigma_xy %*% t(contrast)
        cov2cor(sigma_bxby)[1, 2]
    })
    # keep estimate, lower, upper, p.value
    cor <- cor[, c("estimate", "lower", "upper", "p.value")]
    return(cor)
}

strategy_2b_adjusted <- function(data_longlong, contrast) {
    unstructured_lmm <- lmm(value ~ variable, repetition = ~ variable | id, structure = "UN", data = data_longlong)
    cor <- lava::estimate(unstructured_lmm, function(p){
        sigma_xy <- sigma(unstructured_lmm, p = p)
        sigma_bxby <- contrast %*% sigma_xy %*% t(contrast)
        sigma_bxby_final <- sigma_bxby[1:2, 1:2] - # sigma_slope
                            sigma_bxby[1:2, 3:4] %*%  # sigma_cross
                            solve(sigma_bxby[3:4, 3:4]) %*% # sigma_baseline^(-1)
                            sigma_bxby[3:4, 1:2] # sigma_cross^T
        cov2cor(sigma_bxby_final)[1, 2]
    })
    # keep estimate, lower, upper, p.value
    cor <- cor[, c("estimate", "lower", "upper", "p.value")]
    return(cor)
}

# Uses mmrm to estimate the joint variance covariance matrix. It is faster than LMM, 
# but does not provide confidence intervals for the correlation.
strategy_2b_mmrm <- function(data_longlong, contrast) {
    data_longlong$id <- as.factor(data_longlong$id)
    data_longlong$variable <- as.factor(data_longlong$variable)
    strat2_mmrm <- mmrm::mmrm(value ~ variable + us(variable | id), data = data_longlong)
    sigma_xy <- VarCorr(strat2_mmrm)
    sigma_bxby <- contrast %*% sigma_xy %*% t(contrast)
    cor <- cov2cor(sigma_bxby)[1, 2]
    return(cor)
}

strategy_2b_mmrm_adjusted <- function(data_longlong, contrast) {
    data_longlong$id <- as.factor(data_longlong$id)
    data_longlong$variable <- as.factor(data_longlong$variable)
    strat2_mmrm <- mmrm::mmrm(value ~ variable + us(variable | id), data = data_longlong)
    sigma_xy <- VarCorr(strat2_mmrm)
    sigma_bxby <- contrast %*% sigma_xy %*% t(contrast)
    sigma_bxby_final <- sigma_bxby[1:2, 1:2] - # sigma_slope
                        sigma_bxby[1:2, 3:4] %*%  # sigma_cross
                        solve(sigma_bxby[3:4, 3:4]) %*% # sigma_baseline^(-1)
                        sigma_bxby[3:4, 1:2] # sigma_cross^T

    cor <- cov2cor(sigma_bxby_final)[1, 2]
    return(cor)
}

prepare_data_for_strat_2a <- function(data_longlong) {
    data_long_no_time <- data_longlong %>% dplyr::select(-time)
    data_wide <- dcast(data_long_no_time,
                        id ~ variable,
                        value.var = "value")
    return(data_wide)
}

get_sigma_vb_strat_2a <- function(data, contrast) {
    sigma_vy <- cov(data, use = "pairwise.complete.obs")
    sigma_vb <- contrast %*% sigma_vy %*% t(contrast)
    return(sigma_vb)
}

get_sigma_vb_strat_2a_adjusted <- function(data, contrast) {
    sigma_vy <- cov(data, use = "pairwise.complete.obs")
    sigma_vb <- contrast %*% sigma_vy %*% t(contrast)
    sigma_vb_final <- sigma_vb[1:2, 1:2] - # sigma_slope
                      sigma_vb[1:2, 3:4] %*%  # sigma_cross
                      solve(sigma_vb[3:4, 3:4]) %*% # sigma_baseline^(-1)
                      sigma_vb[3:4, 1:2] # sigma_cross^T
    return(sigma_vb_final)
}

bs_strat_2a <- function(data, indices, contrast) {
    data_boot <- data[indices, ]
    sigma_vb_strat_2a <- get_sigma_vb_strat_2a(data_boot, contrast)
    cor <- cov2cor(sigma_vb_strat_2a)[1, 2]
    return(cor)
}

bs_strat_2a_adjusted <- function(data, indices, contrast) {
    data_boot <- data[indices, ]
    sigma_vb_strat_2a <- get_sigma_vb_strat_2a_adjusted(data_boot, contrast)
    cor <- cov2cor(sigma_vb_strat_2a)[1, 2]
    return(cor)
}

strategy_2a <- function(data, contrast, adjusted=FALSE) {
    data <- copy(data)
    data[, id := NULL]
    if (adjusted) {
        bootstrap <- boot(data = data, stype = "i", statistic = bs_strat_2a_adjusted, R = 1000, contrast = contrast)
    } else {
        bootstrap <- boot(data = data, stype = "i", statistic = bs_strat_2a, R = 1000, contrast = contrast)
    }
    sd_bs <- sd(bootstrap$t, na.rm = TRUE)
    original <- bootstrap$t0
    ci_l <- original - 1.96 * sd_bs
    ci_u <- original + 1.96 * sd_bs
    p.value <- 2*(1-pnorm(abs(original/sd_bs)))
    results <- cbind(original, ci_l, ci_u, p.value)
    colnames(results) <- c("estimate", "lower", "upper", "p.value")
    results <- data.table(results)
    return(results)
}

run_simulation <- function(n_runs, n_subjects, n_measurements, mu, sigma_XY, contrast, adjusted=FALSE, dropout_pattern, mechanism="dropout", beta_mar=0, label="") {
    results_1 <- numeric(n_runs)
    results_2a <- numeric(n_runs)
    results_2b <- numeric(n_runs)

    missingness_fun <- switch(mechanism,
        dropout     = simulate_dropout,
        missingness = simulate_missingness,
        dropout_mar = function(d, p) simulate_dropout_mar(d, p, beta_mar),
        stop("mechanism must be 'dropout', 'missingness' or 'dropout_mar'")
    )

    if (nzchar(label)) message(label)
    pb <- txtProgressBar(min = 0, max = n_runs, style = 3)
    for (i in 1:n_runs) {
        set.seed(i)
        data <- simulate_data(n_subjects, mu, sigma_XY, n_measurements)
        data_long <- make_data_long(data)
        data_missing_long <- missingness_fun(data_long, dropout_pattern)
        data_missing_longlong <- make_data_longlong(data_missing_long)

        if (adjusted) {
            results_1[i] <- strategy_1_adjusted(data_missing_long)$estimate
            results_2a[i] <- strategy_2a(prepare_data_for_strat_2a(data_missing_longlong), contrast, adjusted=adjusted)$estimate
            results_2b[i] <- strategy_2b_mmrm_adjusted(data_missing_longlong, contrast)
        } else {
            results_1[i] <- strategy_1(data_missing_long)$estimate
            results_2a[i] <- strategy_2a(prepare_data_for_strat_2a(data_missing_longlong), contrast)$estimate
            results_2b[i] <- strategy_2b_mmrm(data_missing_longlong, contrast)
        }
        setTxtProgressBar(pb, i)
    }
    close(pb)

    df_results <- data.table(
        variable = rep(c("Strategy 1", "Strategy 2a", "Strategy 2b"), each = n_runs),
        value = c(results_1, results_2a, results_2b)
    )
    return(df_results)
}
