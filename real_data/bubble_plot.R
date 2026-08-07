# =============================================================================
# plot_weights_and_cov.R
#
# Panel 1 & 2 : 5x5 weight matrices (Strategy 1 & 2)
# Panel 3     : 10x10 joint variance/covariance matrix
#
#   - Size   : absolute value
#   - Colour : blue (positive) / red (negative)
#   - Text   : true value printed in centre of bubble
# =============================================================================

library(ggplot2)
library(dplyr)
library(patchwork)

T <- 5

# ── Weight matrices ───────────────────────────────────────────────────────────

W1 <- matrix(c(
  169.07, -56.66, -52.55, -11.85, -14.00,
  -56.66,  65.57,  -0.98,  -3.56,  -4.36,
  -52.55,  -0.98,  50.58,   2.22,   0.73,
  -11.85,  -3.56,   2.22,  11.30,   1.89,
  -14.00,  -4.36,   0.73,   1.89,  15.74
), nrow = T, byrow = TRUE)

W2 <- matrix(c(
   27.16,  13.58,  0.00, -13.58, -27.16,
   13.58,   6.79,  0.00,  -6.79, -13.58,
    0.00,   0.00,  0.00,   0.00,   0.00,
  -13.58,  -6.79,  0.00,   6.79,  13.58,
  -27.16, -13.58,  0.00,  13.58,  27.16
), nrow = T, byrow = TRUE)

# ── Joint variance/covariance matrix (10x10) ──────────────────────────────────

Sigma <- matrix(c(
  14.42, 13.44, 12.78, 12.75, 12.91, -1.61, -1.99, -1.83, -1.58, -1.48,
  13.44, 15.24, 13.33, 13.15, 13.85, -1.15, -1.56, -1.37, -1.26, -0.92,
  12.78, 13.33, 14.85, 13.15, 13.65, -0.22, -0.56, -0.57, -0.13,  0.03,
  12.75, 13.15, 13.15, 15.45, 14.70,  0.02, -0.74, -0.72, -0.67, -0.49,
  12.91, 13.85, 13.65, 14.70, 16.59,  0.04, -0.93, -0.72, -0.40, -0.43,
  -1.61, -1.15, -0.22,  0.02,  0.04, 12.15, 11.08, 11.36, 11.77, 11.86,
  -1.99, -1.56, -0.56, -0.74, -0.93, 11.08, 12.32, 12.01, 12.28, 12.67,
  -1.83, -1.37, -0.57, -0.72, -0.72, 11.36, 12.01, 13.01, 12.64, 13.21,
  -1.58, -1.26, -0.13, -0.67, -0.40, 11.77, 12.28, 12.64, 14.53, 14.05,
  -1.48, -0.92,  0.03, -0.49, -0.43, 11.86, 12.67, 13.21, 14.05, 15.34
), nrow = 2*T, byrow = TRUE)

# ── Labels ────────────────────────────────────────────────────────────────────

labels_T  <- paste0("t[", 0:(T-1), "]")
labels_2T <- c(paste0("X[", 0:(T-1), "]"), paste0("Y[", 0:(T-1), "]"))

# ── Helper: long-format data frame ───────────────────────────────────────────

make_df <- function(M, row_labs, col_labs) {
  n <- nrow(M)
  expand.grid(row = seq_len(n), col = seq_len(n)) |>
    mutate(
      val       = as.vector(M[cbind(row, col)]),
      row_label = factor(row_labs[row], levels = rev(row_labs)),
      col_label = factor(col_labs[col], levels = col_labs)
    )
}

df1    <- make_df(W1,    labels_T,  labels_T)
df2    <- make_df(W2,    labels_T,  labels_T)
df_cov <- make_df(Sigma, labels_2T, labels_2T)

# Size scales: shared for weight panels, separate for covariance
val_max_W   <- max(abs(W1), abs(W2))
val_max_cov <- max(abs(Sigma))

# ── Helper: build plot ────────────────────────────────────────────────────────

make_plot <- function(df,
                      val_max,
                      show_y_labels = TRUE,
                      max_size      = 14,
                      text_size     = 2.2,
                      sep_after     = NULL) {   # draw block separator after this col/row index
  p <- ggplot(df, aes(x = col_label, y = row_label)) +
    geom_point(
      aes(size = abs(val), colour = val > 0),
      shape = 16,
      alpha = 0.45
    ) +
    scale_colour_manual(
      values = c("TRUE" = "#4A90BF", "FALSE" = "#C0392B"),
      guide  = "none"
    ) +
    geom_text(
      aes(label = formatC(val, format = "f", digits = 2)),
      size     = text_size,
      colour   = "grey10",
      fontface = "bold"
    ) +
    scale_size_area(max_size = max_size, limits = c(0, val_max), guide = "none") +
    scale_x_discrete(labels = function(x) parse(text = x), position = "top") +
    scale_y_discrete(
      labels = if (show_y_labels) function(x) parse(text = x)
               else               function(x) rep("", length(x))
    ) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 11) +
    theme(

      axis.text.x      = element_text(size = 9),
      axis.text.y      = element_text(size = 9),
      panel.grid.major = element_line(colour = "grey90"),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(colour = "black", fill = NA),
      legend.position  = "none",
      plot.margin      = margin(10, 6, 10, 6)
    ) +
    coord_fixed()

  # Optional block separator lines (for 2x2 block structure)
  if (!is.null(sep_after)) {
    p <- p +
      geom_vline(xintercept = sep_after + 0.5, colour = "black", linewidth = 0.7) +
      geom_hline(yintercept = nlevels(df$row_label) - sep_after + 0.5,
                 colour = "black", linewidth = 0.7)
  }

  p
}

# ── Build panels ──────────────────────────────────────────────────────────────

p1 <- make_plot(df1,    val_max = val_max_W,   show_y_labels = TRUE)
p2 <- make_plot(df2,    val_max = val_max_W,   show_y_labels = TRUE)
p3 <- make_plot(df_cov, val_max = val_max_cov, show_y_labels = TRUE,
                max_size = 10, text_size = 1.8, sep_after = T)

# ── Layout: weight panels side by side ───────────────────────────────────────

weights <- p1 + p2 + plot_layout(ncol = 2)

print(weights)
print(p3)

# ── Save ──────────────────────────────────────────────────────────────────────

ggsave(
  filename = "weight_strat1_bubble_plot.pdf",
  plot     = p1,
  width    = 6,
  height   = 6,
  units    = "cm"
)

ggsave(
  filename = "weight_strat2_bubble_plot.pdf",
  plot     = p2,
  width    = 6,
  height   = 6,
  units    = "cm"
)

ggsave(
  filename = "cov_bubble_plot.pdf",
  plot     = p3,
  width    = 12,
  height   = 12,
  units    = "cm"
)

message("Plots saved to weight_bubble_plots.pdf and cov_bubble_plot.pdf")