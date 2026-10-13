library(ggplot2)
library(patchwork)
library(data.table)

results <- fread("results/correlation_slope_vs_slope.tsv")
results <- results[failed == FALSE]
results[, error := value - true]
results[, title := factor(paste(scenario, mechanism, sep = " · "),
                          levels = c("S1 · MCAR", "S2 · MCAR", "S3 · MCAR",
                                     "S1b · MAR", "S2b · MAR", "S3b · MAR"))]

# mean +- 2 * MCSE of the error per panel and strategy
summary_error <- results[, .(mean = mean(error), mcse = sd(error) / sqrt(.N)),
                         by = .(adjusted, mechanism, title, strategy)]

strategy_colors <- c("Strategy 1" = "#1b9e77", "Strategy 2a" = "#d95f02", "Strategy 2b" = "#7570b3")

plot_row <- function(adj, mech) {
    ggplot(results[adjusted == adj & mechanism == mech], aes(x = strategy, y = error)) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        geom_boxplot(aes(fill = strategy), alpha = 0.5, outlier.size = 0.5) +
        geom_pointrange(data = summary_error[adjusted == adj & mechanism == mech],
                        aes(y = mean, ymin = mean - 2 * mcse, ymax = mean + 2 * mcse),
                        size = 0.2) +
        facet_wrap(~ title, nrow = 1) +
        scale_fill_manual(name = NULL, values = strategy_colors) +
        theme_minimal() +
        theme(text = element_text(size = 10),
              axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.text.y = element_text(size = 8)) +
        ylab("Estimate - truth")
}

for (adj in c(FALSE, TRUE)) {
    p_all <- (plot_row(adj, "MCAR") / plot_row(adj, "MAR")) +
        plot_layout(guides = "collect") & theme(legend.position = "bottom")
    print(p_all)

    out_file <- paste0("results/simulation_bias_", ifelse(adj, "adjusted", "unadjusted"))
    ggsave(paste0(out_file, ".pdf"), p_all, width = 16, height = 12, units = "cm")
    ggsave(paste0(out_file, ".png"), p_all, width = 16, height = 12, dpi = 300, units = "cm", bg = "white")
}
