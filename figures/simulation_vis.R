library(ggplot2)

results <- read.table("results/correlation_slope_vs_slope.tsv", sep = "\t", header = TRUE)

results$adjusted_text <- ifelse(results$adjusted, "Adjusted", "Unadjusted")
results$experiment_adjusted <- paste(results$experiment, results$adjusted_text, sep = "\n")
results$experiment_adjusted <- factor(results$experiment_adjusted,
                                     levels = c("Scenario 1\nUnadjusted", "Scenario 2\nUnadjusted", "Scenario 3\nUnadjusted",
                                               "Scenario 1\nAdjusted", "Scenario 2\nAdjusted", "Scenario 3\nAdjusted"))

p_all <- ggplot(results, aes(x = strategy, y = value)) +
    geom_boxplot() +
    facet_wrap(~ experiment_adjusted, scales = "free_y", ncol = 3) +
    theme_minimal() +
    geom_hline(aes(yintercept = true, linetype = "True", color = "True")) +
    scale_linetype_manual(name = NULL, values = c("True" = "dotted")) +
    scale_color_manual(name = NULL, values = c("True" = "blue")) +
    guides(linetype = guide_legend(title = NULL), color = guide_legend(title = NULL)) +
    theme(text = element_text(size = 10)) +
    ylab("Correlation")

p_all <- p_all + 
    theme(axis.text.x = element_text(size = 8, hjust = 1, angle = 45), 
          axis.title.x = element_blank(), 
          axis.text.y = element_text(size = 8),
          legend.position = "none")
print(p_all)

ggsave("results/correlation_slope_vs_slope.pdf", p_all, width = 12, height = 9, dpi = 300, units = "cm")
