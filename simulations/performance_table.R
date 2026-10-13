library(data.table)

# Performance table per scenario x strategy, computed from the simulation results.
# n_sim is the number of successful runs (failed runs are excluded from all metrics).
results <- fread("results/correlation_slope_vs_slope.tsv")

scenario_levels <- c("S1", "S2", "S3", "S1b", "S2b", "S3b")

fmt <- function(x) sprintf("%.3f", x)

for (adj in c(FALSE, TRUE)) {
    performance <- results[adjusted == adj, {
        est   <- value[!failed]
        n_sim <- length(est)
        emp_se <- sd(est)
        list(
            true_rho  = fmt(true[1]),
            mean_rho  = fmt(mean(est)),
            bias      = paste0(fmt(mean(est - true[1])), " (", fmt(emp_se / sqrt(n_sim)), ")"),
            emp_se    = paste0(fmt(emp_se), " (", fmt(emp_se / sqrt(2 * (n_sim - 1))), ")"),
            rmse      = fmt(sqrt(mean((est - true[1])^2))),
            failures  = sum(failed)
        )
    }, by = .(scenario, mechanism, strategy)]
    performance <- performance[order(factor(scenario, levels = scenario_levels), strategy)]
    print(performance)

    out_file <- paste0("results/performance_", ifelse(adj, "adjusted", "unadjusted"), ".tsv")
    write.table(performance, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
}
