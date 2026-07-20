source("R/functions.R")

library(readxl)
library(dplyr)
library(data.table)
library(ggplot2)
library(patchwork)
library(LMMstar)
library(prodlim)
library(riskRegression)
library(lme4)
library(lmerTest)

set.seed(42)

# NOTE: data/data.xlsx is not publicly available.
# This script documents the analysis applied to the real ADNI dataset.

data <- read_excel("data/data.xlsx")
print(dim(data))

data_keep <- data[data$keep == 1, ]
print(dim(data_keep))

data_keep <- data_keep %>%
    dplyr::mutate(
        subject_id = factor(subject_id),
        session_id = factor(session_id),
        mri_field  = factor(mri_field),
        gender     = factor(gender),
        diagnosis  = factor(diagnosis, levels = c("CN", "MCI", "AD"))
    )

data_keep <- data.table(data_keep)

data_subset <- data_keep %>% dplyr::select(
    subject_id,
    session_id,
    diagnosis,
    pad_pyment,
    gm_icv,
    ADNI_MEM
)

data_subset <- data_subset %>%
    dplyr::mutate(
        subject_id = factor(subject_id),
        session_id = factor(session_id),
        diagnosis  = factor(diagnosis, levels = c("CN", "MCI", "AD"))
    )

data_subset$time <- as.numeric(gsub("M", "", data_subset$session_id, fixed = TRUE))
data_subset$time <- data_subset$time / 12

data_4years <- data_subset %>% dplyr::filter(time <= 4)
print("Number of scans until 4 years:")
print(dim(data_4years))

data_4years_min2 <- data_4years %>%
    dplyr::group_by(subject_id) %>%
    dplyr::filter(n() >= 2) %>%
    dplyr::ungroup()
print("Number of subjects with at least 2 visits:")
print(dim(data_4years_min2))

print("Scans per session_id:")
data_4years_min2 <- droplevels(data_4years_min2)
print(table(data_4years_min2$session_id))

# Only keep yearly visits (drop 6-/18-/30-/42-month interim scans)
data_4years_min2_yearly <- data_4years_min2[data_4years_min2$session_id %in% c("M006", "M018", "M030", "M042") == FALSE, ]
data_4years_min2_yearly <- droplevels(data_4years_min2_yearly)
print(table(data_4years_min2_yearly$session_id))

table_data <- table(data_4years_min2_yearly$session_id)
print("Drop-out percentage:")
print((1 - table_data / table_data[1]) * 100)

# Restrict to CN-at-baseline subjects
diag    <- "CN"
data_cn <- data_4years_min2_yearly %>%
    dplyr::arrange(subject_id, time) %>%
    dplyr::group_by(subject_id) %>%
    dplyr::filter(dplyr::first(diagnosis) == diag) %>%
    dplyr::ungroup()

# Normalise time so that each subject's first visit is t = 0
data_cn <- data_cn %>%
    dplyr::group_by(subject_id) %>%
    dplyr::mutate(time = time - min(time)) %>%
    dplyr::ungroup()

print(dim(data_cn))
print(table(data_cn$session_id))

data_ready <- data_cn %>%
    dplyr::select(subject_id, time, pad_pyment, gm_icv) %>%
    dplyr::rename(id = subject_id, X = pad_pyment, Y = gm_icv)

# ---------------------------------------------------------------------------
# Trajectories
# ---------------------------------------------------------------------------

p_diag_X <- ggplot(data_ready, aes(x = time, y = X, group = id)) +
    geom_line(alpha = 0.1, color = "#3b3b3b") +
    geom_point(color = "#0f0f0f", size = 0.5, alpha = 0.5) +
    labs(x = "Time", y = "PAD") +
    theme_minimal() +
    theme(axis.text.x  = element_text(size = 8),
          axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          axis.text.y  = element_text(size = 8),
          title          = element_text(size = 10),
          legend.title   = element_blank(),
          legend.position = "none")

p_diag_Y <- ggplot(data_ready, aes(x = time, y = Y, group = id)) +
    geom_line(alpha = 0.1, color = "#3b3b3b") +
    geom_point(color = "#0f0f0f", size = 0.5, alpha = 0.5) +
    labs(x = "Time", y = "Gray Matter Volume") +
    theme_minimal() +
    theme(axis.text.x  = element_text(size = 8),
          axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          axis.text.y  = element_text(size = 8),
          title          = element_text(size = 10),
          legend.title   = element_blank(),
          legend.position = "none")

p_diag <- p_diag_X + p_diag_Y + plot_layout(nrow = 1)
ggsave("results/data_real.pdf", p_diag, width = 12, height = 12/2.6, dpi = 300, units = "cm")

# ---------------------------------------------------------------------------
# Descriptives: mean per timepoint
# ---------------------------------------------------------------------------

print("Mean at each timepoint:")
print(data_ready %>%
    dplyr::group_by(time) %>%
    dplyr::summarise(
        n        = n(),
        mean_gm  = mean(Y, na.rm = TRUE),
        mean_pad = mean(X, na.rm = TRUE)
    ))

# ---------------------------------------------------------------------------
# Joint covariance matrix
# ---------------------------------------------------------------------------

data_ready   <- as.data.table(data_ready)
data_longlong <- make_data_longlong(data_ready)

lmm_covvar <- lmm(
    value ~ variable,
    repetition = ~ variable | id,
    structure  = "UN",
    data       = data_longlong
)
sigma_lmm <- sigma(lmm_covvar)
print("Estimated joint covariance matrix:")
print(round(sigma_lmm, 2))

# ---------------------------------------------------------------------------
# Correlation estimation: all three strategies
# ---------------------------------------------------------------------------

n_measurements <- length(unique(data_ready$time))
time    <- unique(data_ready$time)
weight  <- (time - mean(time)) / sum((time - mean(time))^2)
contrast <- rbind(c(weight, rep(0, n_measurements)), c(rep(0, n_measurements), weight))

contrast_adjusted <- rbind(
    c(weight, rep(0, n_measurements)),
    c(rep(0, n_measurements), weight),
    c(1, rep(0, n_measurements - 1), rep(0, n_measurements)),
    c(rep(0, n_measurements), c(1, rep(0, n_measurements - 1)))
)

results_1          <- strategy_1(data_ready)
results_2a         <- strategy_2a(prepare_data_for_strat_2a(data_longlong), contrast)
results_2b         <- strategy_2b(data_longlong, contrast)
results_1_adjusted  <- strategy_1_adjusted(data_ready)
results_2a_adjusted <- strategy_2a(prepare_data_for_strat_2a(data_longlong), contrast_adjusted, adjusted = TRUE)
results_2b_adjusted <- strategy_2b_adjusted(data_longlong, contrast_adjusted)

results_all <- rbind(
    cbind(strategy = "1",           results_1),
    cbind(strategy = "2a",          results_2a),
    cbind(strategy = "2b",          results_2b),
    cbind(strategy = "1 adjusted",  results_1_adjusted),
    cbind(strategy = "2a adjusted", results_2a_adjusted),
    cbind(strategy = "2b adjusted", results_2b_adjusted)
)

results_all <- results_all %>%
    dplyr::mutate(
        estimate = round(estimate, 4),
        lower    = round(lower,    4),
        upper    = round(upper,    4),
        p.value  = round(p.value,  4)
    )

print("All results:")
print(results_all)
write.csv(results_all, "results/real_data_results.csv", row.names = FALSE)

# ---------------------------------------------------------------------------
# Cross-timepoint weight product tables
#
# Left  table : sum_k  n_k w_i^k * w_j^k     for each pair (i, j)
# Right table : normalised by sum_i sum_j sum_k |n_k w_i^k * w_j^k|
# ---------------------------------------------------------------------------

n_ready <- nrow(data_ready[data_ready$time == 0, ])

time_patterns <- data_ready %>%
    group_by(id) %>%
    summarise(time_combo = list(sort(unique(time)))) %>%
    group_by(time_combo) %>%
    summarise(
        n_subjects = n(),
        times      = paste(time_combo[[1]], collapse = ", "),
        weights    = list({
            t <- time_combo[[1]]
            if (length(t) == 1) {
                1  # weight = 1 for single time points
            } else {
                (t - mean(t)) / sum((t - mean(t))^2)
            }
        }),
        weighted_weights = list(weights[[1]] * n_subjects)
    ) %>%
    arrange(desc(n_subjects)) %>%
    select(n_subjects, times, weights, weighted_weights)

all_times <- sort(unique(unlist(strsplit(paste(time_patterns$times, collapse = ", "), ", "))))
all_times <- as.numeric(all_times)

K <- nrow(time_patterns)
T <- length(all_times)

cross_weight_matrix <- matrix(0, nrow = T, ncol = T, dimnames = list(all_times, all_times))

for (k in 1:K) {
    pattern_times <- as.numeric(strsplit(time_patterns$times[k], ", ")[[1]])
    idx <- match(pattern_times, all_times)
    wk  <- time_patterns$weights[[k]]
    nk  <- time_patterns$n_subjects[k]
    cross_weight_matrix[idx, idx] <- cross_weight_matrix[idx, idx] + nk * outer(wk, wk)
}

cat("\nLeft table: sum_k n_k * w_i^k * w_j^k\n")
print(round(cross_weight_matrix, 2))

abs_sum <- 0
for (k in 1:K) {
    wk      <- time_patterns$weights[[k]]
    nk      <- time_patterns$n_subjects[k]
    abs_sum <- abs_sum + nk * sum(abs(outer(wk, wk)))
}

normalised_matrix <- cross_weight_matrix / abs_sum
cat("\nRight table: normalised sum_k n_k * w_i^k * w_j^k\n")
print(round(normalised_matrix * 100, 2))

# ---------------------------------------------------------------------------
# Sanity check
#
# Tile cross_weight_matrix into a 2x2 block, elementwise-multiply with
# sigma_lmm, then sum each block:
#   top-left  -> var(beta_x),          top-right -> cov(beta_x, beta_y)
#   bot-right -> var(beta_y)
# The implied correlation should match strategy 2b.
# ---------------------------------------------------------------------------

W            <- cross_weight_matrix
block_matrix <- rbind(cbind(W, W), cbind(W, W))
ep           <- block_matrix * sigma_lmm

var_x  <- sum(ep[1:T,         1:T        ])
cov_xy <- sum(ep[1:T,        (T+1):(2*T) ])
var_y  <- sum(ep[(T+1):(2*T), (T+1):(2*T)])

corr_check <- cov_xy / sqrt(var_x * var_y)
cat("\nSanity check (observed dropout pattern):\n")
cat("  var(beta_x)         :", round(var_x,       6), "\n")
cat("  cov(beta_x, beta_y) :", round(cov_xy,      6), "\n")
cat("  var(beta_y)         :", round(var_y,        6), "\n")
cat("  implied correlation :", round(corr_check,   6), "\n")
cat("  strategy 1 estimate :", round(results_1$estimate, 6), "\n")

# ---------------------------------------------------------------------------
# Ideal weight matrix (no missing data)
#
# Under complete data every subject uses the same weights, so
# W_ideal[i,j] = n_ready * weight[i] * weight[j].
# Repeating the sanity check with W_ideal should recover strategy 2b.
# ---------------------------------------------------------------------------

W_ideal           <- n_ready * outer(weight, weight)
rownames(W_ideal) <- all_times
colnames(W_ideal) <- all_times

cat("\nIdeal cross_weight_matrix (complete data):\n")
print(round(W_ideal, 2))

block_ideal <- rbind(cbind(W_ideal, W_ideal), cbind(W_ideal, W_ideal))
ep_ideal    <- block_ideal * sigma_lmm

var_x_ideal  <- sum(ep_ideal[1:T,         1:T        ])
cov_xy_ideal <- sum(ep_ideal[1:T,        (T+1):(2*T) ])
var_y_ideal  <- sum(ep_ideal[(T+1):(2*T), (T+1):(2*T)])

corr_ideal <- cov_xy_ideal / sqrt(var_x_ideal * var_y_ideal)
cat("\nSanity check (ideal / complete data):\n")
cat("  var(beta_x)         :", round(var_x_ideal,  6), "\n")
cat("  cov(beta_x, beta_y) :", round(cov_xy_ideal, 6), "\n")
cat("  var(beta_y)         :", round(var_y_ideal,  6), "\n")
cat("  implied correlation :", round(corr_ideal,   6), "\n")
cat("  strategy 2b estimate:", round(results_2b$estimate, 6), "\n")

# ---------------------------------------------------------------------------
# Survival analysis: probability of remaining CN in the cohort
# ---------------------------------------------------------------------------

data_surv <- data_cn %>%
    dplyr::group_by(subject_id) %>%
    dplyr::arrange(time) %>%
    dplyr::summarise(
        time4 = ifelse(
            any(diagnosis != "CN"),
            min(time[diagnosis != "CN"]),
            max(time)
        ),
        diag4 = as.integer(any(diagnosis != "CN"))
    ) %>%
    dplyr::ungroup()

print("Conversion summary (0 = censored, 1 = converted from CN):")
print(table(data_surv$diag4))

e.KM <- prodlim(Hist(time4, diag4) ~ 1, data = data_surv)

pdf("results/surv_CN.pdf", width = 2.5, height = 3.5)
print(
    ggprodlim(e.KM, position_atrisk = FALSE) +
        ggtitle("CN cohort") +
        labs(y = "Probability of remaining CN")
)
dev.off()

# ---------------------------------------------------------------------------
# Linear mixed models: PAD ~ GM volume
# ---------------------------------------------------------------------------

lmm_ri <- lmer(X ~ Y + (1 | id),    data = data_ready)
lmm_rs <- lmer(X ~ Y + (time | id), data = data_ready)

print("Random intercept model: PAD ~ volume + (1|id)")
print(summary(lmm_ri))

print("Random slope model: PAD ~ volume + (time|id)")
print(summary(lmm_rs))
