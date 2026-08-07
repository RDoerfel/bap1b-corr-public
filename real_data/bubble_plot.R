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
  166.24, -56.82, -51.27, -11.15, -14.00,
  -56.82,  66.35,  -0.82,  -3.66,  -5.06,
  -51.27,  -0.82,  50.36,   2.18,  -0.45,
  -11.15,  -3.66,   2.18,  10.71,   1.92,
  -14.00,  -5.06,  -0.45,   1.92,  18.60
), nrow = T, byrow = TRUE)

W2 <- matrix(c(
  26.64,  13.32, 0, -13.32, -26.64,
  13.32,   6.66, 0,  -6.66, -13.32,
  0.00,   0.00, 0,   0.00,   0.00,
  -13.32,  -6.66, 0,   6.66,  13.32,
  -26.64, -13.32, 0,  13.32,  26.64
), nrow = T, byrow = TRUE)

# ── Joint variance/covariance matrix (10x10) ──────────────────────────────────

Sigma <- matrix(c(
  14.45, 13.39, 12.83, 12.75, 12.88, -1.64, -1.81, -1.88, -1.62, -1.40,
  13.39, 15.14, 13.32, 13.09, 13.84, -1.08, -1.33, -1.37, -1.24, -0.81,
  12.83, 13.32, 14.80, 13.14, 13.62, -0.34, -0.55, -0.74, -0.36, -0.04,
  12.75, 13.09, 13.14, 15.32, 14.48, -0.50, -0.49, -0.97, -0.89, -0.48,
  12.88, 13.84, 13.62, 14.48, 16.69,  0.08, -0.44, -0.72, -0.21, -0.21,
  -1.64, -1.08, -0.34, -0.50,  0.08, 12.08, 11.15, 11.35, 11.80, 11.82,
  -1.81, -1.33, -0.55, -0.49, -0.44, 11.15, 12.35, 11.89, 12.30, 12.58,
  -1.88, -1.37, -0.74, -0.97, -0.72, 11.35, 11.89, 13.04, 12.61, 13.15,
  -1.62, -1.24, -0.36, -0.89, -0.21, 11.80, 12.30, 12.61, 14.63, 14.02,
  -1.40, -0.81, -0.04, -0.48, -0.21, 11.82, 12.58, 13.15, 14.02, 15.20
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
  filename = "results/weight_strat1_bubble_plot.pdf",
  plot     = p1,
  width    = 6,
  height   = 6,
  units    = "cm"
)

ggsave(
  filename = "results/weight_strat2_bubble_plot.pdf",
  plot     = p2,
  width    = 6,
  height   = 6,
  units    = "cm"
)

ggsave(
  filename = "results/cov_bubble_plot.pdf",
  plot     = p3,
  width    = 12,
  height   = 12,
  units    = "cm"
)

message("Plots saved to weight_bubble_plots.pdf and cov_bubble_plot.pdf")