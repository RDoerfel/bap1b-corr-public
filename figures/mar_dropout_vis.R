library(ggplot2)

results <- read.table("results/mar_dropout.tsv", sep = "\t", header = TRUE)

results$experiment <- factor(results$experiment,
                             levels = c("MCAR", "MAR beta=0", "MAR beta=-1", "MAR beta=-2"),
                             labels = c("MCAR dropout", "MAR dropout\nbeta = 0", "MAR dropout\nbeta = -1", "MAR dropout\nbeta = -2"))

p_all <- ggplot(results, aes(x = strategy, y = value)) +
    geom_boxplot() +
    facet_wrap(~ experiment, ncol = 4) +
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

ggsave("results/mar_dropout.pdf", p_all, width = 16, height = 7, dpi = 300, units = "cm")
