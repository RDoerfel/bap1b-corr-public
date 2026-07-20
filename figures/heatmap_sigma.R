# Define joint variance-covariance structures for the three simulation scenarios.
n_measurements <- 5
sigma_marginalX <- rep(1, n_measurements)
sigma_marginalY <- rep(1, n_measurements)
sigma_marginalXY <- diag(c(sigma_marginalX, sigma_marginalY))

# Scenario 1: diagonal cross-correlation
cor_X <- 1
cor_Y <- 1
cor_XY_diag <- 0.6

rho_X <- diag(cor_X, n_measurements, n_measurements)
rho_Y <- diag(cor_Y, n_measurements, n_measurements)
rho_XY_diagonal <- diag(cor_XY_diag, n_measurements, n_measurements)

# combine the full correlation matrix
rho_diag <- rbind(cbind(rho_X, t(rho_XY_diagonal)), cbind(rho_XY_diagonal, rho_Y))
sigma_XY_diag_const <- sigma_marginalXY %*% rho_diag %*% sigma_marginalXY

# Scenario 2: diagonal cross-correlation non-constant cross-correlation
cor_X <- 1
cor_Y <- 1
cor_XY_diag <- c(0.8, 0.7, 0.6, 0.5, 0.4)

rho_X <- diag(cor_X, n_measurements, n_measurements)
rho_Y <- diag(cor_Y, n_measurements, n_measurements)
rho_XY_diagonal <- diag(cor_XY_diag)

# combine the full correlation matrix
rho_diag <- rbind(cbind(rho_X, t(rho_XY_diagonal)), cbind(rho_XY_diagonal, rho_Y))
sigma_XY_diag <- sigma_marginalXY %*% rho_diag %*% sigma_marginalXY

# Scenario 3: realistic covariance structure, estimated from real data (ADNI, longitudinal, CN)
rho_YY <- matrix(
    c(
        1.00, 0.91, 0.89, 0.87, 0.86,
        0.91, 1.00, 0.94, 0.92, 0.92,
        0.89, 0.94, 1.00, 0.91, 0.95,
        0.87, 0.92, 0.91, 1.00, 0.94,
        0.86, 0.92, 0.95, 0.94, 1.00
    ),
    nrow = 5, ncol = 5, byrow = TRUE,
    dimnames = list(
        paste0("Y_", c("0", "1", "2", "3", "4")),
        paste0("Y_", c("0", "1", "2", "3", "4"))
    )
)

# X (pad) correlation matrix
rho_XX <- matrix(
    c(
        1.00, 0.91, 0.88, 0.86, 0.81,
        0.91, 1.00, 0.89, 0.86, 0.86,
        0.88, 0.89, 1.00, 0.88, 0.89,
        0.86, 0.86, 0.88, 1.00, 0.91,
        0.81, 0.86, 0.89, 0.91, 1.00
    ),
    nrow = 5, ncol = 5, byrow = TRUE,
    dimnames = list(
        paste0("X_", c("0", "1", "2", "3", "4")),
        paste0("X_", c("0", "1", "2", "3", "4"))
    )
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
    dimnames = list(
        paste0("Y_", c("0", "1", "2", "3", "4")),
        paste0("X_", c("0", "1", "2", "3", "4"))
    )
)

# Cross-correlation block: X (pad) rows, Y (gm) columns (transpose of above)
rho_XY <- t(rho_YX)
rho <- rbind(
    cbind(rho_XX, rho_XY),
    cbind(t(rho_XY), rho_YY)
)

sigma_marginalXY <- diag(c(rep(1, n_measurements), rep(1, n_measurements)))
sigma_joint_XY <- sigma_marginalXY %*% rho %*% sigma_marginalXY

library(ggplot2)
library(reshape2)
library(gridExtra)

# Function to create heatmap for a covariance matrix (no legend)

plot_cov_heatmap <- function(cov_mat, title = "Covariance Matrix") {
  melted_cov <- melt(cov_mat)
  
  # Create labels - same for both axes
  n <- nrow(cov_mat)
  labels <- c("X1", "X2", "X3", "X4", "X5", "Y1", "Y2", "Y3", "Y4", "Y5")[1:n]
  
  # Convert to factors with custom labels
  melted_cov$Var1 <- factor(melted_cov$Var1, levels = 1:n, labels = labels)
  melted_cov$Var2 <- factor(melted_cov$Var2, levels = rev(1:n), labels = rev(labels))
  
  # Create custom color mapping
  colors <- colorRampPalette(c("#2C7BB6", "#ABD9E9", "white", "#FDAE61", "#D7191C"))(100)
  
  ggplot(melted_cov, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradientn(colors = colors, limits = c(-1, 1), name = NULL,
                         guide = guide_colorbar(barwidth = 0.5, barheight = 10)) +
    coord_fixed() +
    theme_minimal() +
    labs(title = title, x = "", y = "") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5))
}


# Create plots
p1 <- plot_cov_heatmap(sigma_XY_diag_const, "Scenario 1")
p2 <- plot_cov_heatmap(sigma_XY_diag, "Scenario 2")
p3 <- plot_cov_heatmap(sigma_joint_XY, "Scenario 3")

# Extract legend from p3
get_legend <- function(plot) {
  tmp <- ggplot_gtable(ggplot_build(plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  tmp$grobs[[leg]]
}

legend <- get_legend(p3)

# Remove legends from all plots
p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")
p3 <- p3 + theme(legend.position = "none")

combined_plot <- grid.arrange(p1, p2, p3, legend, ncol = 4, widths = c(1, 1, 1, 0.3))
ggsave("results/heatmap_covariance_matrices.pdf", combined_plot, width = 17, height = 7, unit="cm")
