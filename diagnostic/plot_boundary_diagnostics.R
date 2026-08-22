args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("usage: Rscript diagnostic/plot_boundary_diagnostics.R diagnostics.csv")
}

results <- read.csv(args[[1]])
results$pdf_ratio <- results$estimated_pdf / results$true_pdf
results$pdf_ratio[!is.finite(results$pdf_ratio) | results$pdf_ratio <= 0] <-
  NA_real_
results$cdf_error <- abs(results$estimated_cdf - results$true_cdf)

tail_results <- subset(results, position <= 1e-2 & true_pdf > 0)
summary <- aggregate(
  cbind(log_pdf_error = abs(log(pdf_ratio)), cdf_error = cdf_error) ~
    scenario + scale + degree + grid_size + bandwidth,
  data = tail_results,
  FUN = function(x) max(x, na.rm = TRUE)
)
summary <- summary[order(summary$log_pdf_error, summary$cdf_error), ]
cat("Worst tail errors\n")
print(head(summary[order(-summary$log_pdf_error), ], 20L), row.names = FALSE)

results$scaled_pdf <- results$estimated_pdf * results$scale
scale_summary <- aggregate(
  cbind(scaled_pdf = scaled_pdf, estimated_cdf = estimated_cdf) ~
    scenario + degree + grid_size + bandwidth + position,
  data = subset(results, true_pdf > 0),
  FUN = function(x) max(x) - min(x)
)
scale_summary$discrepancy <- pmax(
  scale_summary$scaled_pdf,
  scale_summary$estimated_cdf
)
cat("\nLargest scale-equivariance discrepancies\n")
print(
  head(scale_summary[order(-scale_summary$discrepancy), ], 20L),
  row.names = FALSE
)

pdf(sub("[.]csv$", ".pdf", args[[1]]))
for (scenario_name in unique(results$scenario)) {
  selected <- subset(
    results,
    scenario == scenario_name & scale == 1 & grid_size == 400 &
      bandwidth == "0.300000" & position > 0 & position <= 0.1
  )
  matplot(
    selected$position[selected$degree == 0],
    sapply(0:2, function(polynomial_degree) {
      selected$estimated_pdf[selected$degree == polynomial_degree]
    }),
    type = "l",
    log = "x",
    lty = 1,
    xlab = "distance from boundary / support scale",
    ylab = "density",
    main = scenario_name
  )
  lines(
    selected$position[selected$degree == 0],
    selected$true_pdf[selected$degree == 0],
    lty = 2,
    lwd = 2
  )
  legend("topright", c("degree 0", "degree 1", "degree 2", "truth"),
         col = c(1:3, 1), lty = c(1, 1, 1, 2), bty = "n")
}
dev.off()
