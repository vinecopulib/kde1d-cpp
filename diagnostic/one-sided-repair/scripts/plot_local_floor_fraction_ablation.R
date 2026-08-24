# Plot the local bandwidth-floor fraction ablation.
# Run from the root of the kde1d R package.

if (!requireNamespace("pkgload", quietly = TRUE) || !file.exists("DESCRIPTION")) {
  stop("run from the kde1d package root with pkgload installed")
}
pkgload::load_all(".", quiet = TRUE)
source(paste0(
  "inst/include/kde1d-cpp/diagnostic/one-sided-repair/scripts/",
  "one_sided_estimator.R"
))

input_file <- paste0(
  "inst/include/kde1d-cpp/diagnostic/one-sided-repair/results/",
  "local-floor-fraction-ablation.csv"
)
output_file <- paste0(
  "inst/include/kde1d-cpp/diagnostic/one-sided-repair/plots/",
  "local-floor-fraction-ablation.pdf"
)
results <- read.csv(input_file)
fractions <- as.numeric(sub(
  "pct$", "", sub("fixed-", "", results$method[grepl("^fixed", results$method)]),
  fixed = FALSE
))
fractions <- sort(unique(fractions))
comparison_distributions <- list(
  exponential = list(random = rexp, density = dexp, quantile = qexp),
  gamma_2 = list(
    random = function(sample_size) rgamma(sample_size, 2),
    density = function(values) dgamma(values, 2),
    quantile = function(probabilities) qgamma(probabilities, 2)
  ),
  lomax_2 = list(
    random = function(sample_size) (1 - runif(sample_size))^(-0.5) - 1,
    density = function(values) 2 / (1 + values)^3,
    quantile = function(probabilities) (1 - probabilities)^(-0.5) - 1
  )
)
evaluation_probabilities <- sort(unique(c(
  10^seq(-10, -2, length.out = 101L),
  seq(0.01, 0.999, length.out = 401L),
  0.9999, 0.99999
)))
method_colors <- c(
  current = "grey45", `25%` = "#0072B2", `75%` = "#009E73",
  `100%` = "#D55E00"
)

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
pdf(output_file, width = 11, height = 7.5)
par(mfrow = c(2, 2), mar = c(4.3, 4.5, 2.3, 0.8), oma = c(0, 0, 3, 0))
for (sample_size in sort(unique(results$sample_size))) {
  current_results <- results[
    results$sample_size == sample_size & results$method == "current",
  ]
  fixed_results <- results[
    results$sample_size == sample_size & grepl("^fixed", results$method),
  ]
  fixed_results <- merge(
    fixed_results,
    current_results[, c("distribution", "global_ise", "boundary_ise")],
    by = "distribution", suffixes = c("", "_current")
  )
  for (metric in c("global", "boundary")) {
    ratio_matrix <- vapply(fractions, function(fraction) {
      selected_results <- fixed_results[
        fixed_results$method == paste0("fixed-", fraction, "pct"),
      ]
      selected_results[[paste0(metric, "_ise")]] /
        selected_results[[paste0(metric, "_ise_current")]]
    }, numeric(nrow(current_results)))
    matplot(
      fractions, t(ratio_matrix), type = "l", lty = 1, lwd = 1.4,
      col = adjustcolor("grey40", 0.45), ylim = range(c(0.7, 1.1, ratio_matrix)),
      xlab = "lowest sample fraction used for floor (%)",
      ylab = "candidate / current ISE",
      main = sprintf("%s ISE, n=%d", metric, sample_size)
    )
    lines(fractions, colMeans(ratio_matrix), lwd = 3, col = "#009E73")
    abline(h = 1, lty = 2)
  }
}
mtext("Bandwidth-floor fraction: density-specific curves and mean", outer = TRUE,
      font = 2, cex = 1.1)

for (boundary_zoom in c(FALSE, TRUE)) {
  par(mfrow = c(2, 3), mar = c(4.2, 4.2, 2.5, 0.5), oma = c(0, 0, 3, 0))
  for (sample_size in c(100L, 1000L)) {
    for (distribution_index in seq_along(comparison_distributions)) {
      evaluation_points <- comparison_distributions[[distribution_index]]$quantile(
        evaluation_probabilities
      )
      true_density <- comparison_distributions[[distribution_index]]$density(
        evaluation_points
      )
      set.seed(
        20260914L + 100000L * match(
          names(comparison_distributions)[[distribution_index]],
          c("exponential", "half_normal", "gamma_2", "gamma_0.75", "lomax_2")
        ) + 1000L * match(sample_size, c(100L, 1000L)) + 1L
      )
      observations <- comparison_distributions[[distribution_index]]$random(
        sample_size
      )
      estimates <- list(
        current = normalize_one_sided_density(
          evaluation_points,
          dkde1d(evaluation_points, kde1d(observations, xmin = 0))
        ),
        `25%` = fit_one_sided_candidate(
          observations, evaluation_points, local_floor_fraction = 0.25
        )$density,
        `75%` = fit_one_sided_candidate(
          observations, evaluation_points, local_floor_fraction = 0.75
        )$density,
        `100%` = fit_one_sided_candidate(
          observations, evaluation_points, local_floor_fraction = 1
        )$density
      )
      plotted_points <- if (boundary_zoom) {
        evaluation_probabilities <= 0.2
      } else {
        evaluation_probabilities <= 0.999
      }
      plot(
        evaluation_points[plotted_points], true_density[plotted_points],
        type = "l", lwd = 3, col = "black", log = if (boundary_zoom) "x" else "",
        ylim = c(0, quantile(c(
          true_density[plotted_points],
          unlist(lapply(estimates, function(density) density[plotted_points]))
        ), 0.995)),
        xlab = if (boundary_zoom) "distance from boundary" else "x",
        ylab = "density",
        main = sprintf(
          "%s, n=%d", names(comparison_distributions)[[distribution_index]],
          sample_size
        )
      )
      for (method in names(estimates)) {
        lines(
          evaluation_points[plotted_points], estimates[[method]][plotted_points],
          col = method_colors[[method]], lwd = if (method == "75%") 2.5 else 1.5,
          lty = match(method, names(estimates))
        )
      }
      if (sample_size == 100L && distribution_index == 1L) {
        legend(
          "topright", c("truth", names(estimates)),
          col = c("black", method_colors), lty = c(1, seq_along(estimates)),
          lwd = c(3, 1.5, 1.5, 2.5, 1.5), bty = "n", cex = 0.72
        )
      }
    }
  }
  mtext(
    if (boundary_zoom) "Representative boundary views" else
      "Representative full-density views",
    outer = TRUE, font = 2, cex = 1.1
  )
}
dev.off()
message("wrote ", normalizePath(output_file))
