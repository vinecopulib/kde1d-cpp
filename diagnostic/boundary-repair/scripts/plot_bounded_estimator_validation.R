# Plot the paired validation of the selected bounded estimator.
# Run from the root of the kde1d R package.

if (!requireNamespace("pkgload", quietly = TRUE) || !file.exists("DESCRIPTION")) {
  stop("run from the kde1d package root with pkgload installed")
}
pkgload::load_all(".", quiet = TRUE)
source(paste0(
  "inst/include/kde1d-cpp/diagnostic/boundary-repair/scripts/",
  "bounded_estimator.R"
))

arguments <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(arguments) >= 1L) arguments[[1L]] else paste0(
  "inst/include/kde1d-cpp/diagnostic/boundary-repair/results/",
  "bounded-estimator-validation.csv"
)
output_file <- if (length(arguments) >= 2L) arguments[[2L]] else paste0(
  "inst/include/kde1d-cpp/diagnostic/boundary-repair/plots/",
  "bounded-estimator-validation.pdf"
)
if (length(arguments) > 2L) {
  stop("usage: plot_bounded_estimator_validation.R [input.csv] [output.pdf]")
}

results <- read.csv(input_file)
replication_count <- unique(results$replications)
if (length(replication_count) != 1L) {
  stop("input must contain a single replication count")
}
method_colors <- c(current = "grey45", selected = "#009E73")
validation_distributions <- list(
  uniform = list(random = runif, density = dunif),
  `beta-1-2` = list(
    random = function(sample_size) rbeta(sample_size, 1, 2),
    density = function(values) dbeta(values, 1, 2)
  ),
  `beta-0.75-2` = list(
    random = function(sample_size) rbeta(sample_size, 0.75, 2),
    density = function(values) dbeta(values, 0.75, 2)
  ),
  `beta-0.75-0.75` = list(
    random = function(sample_size) rbeta(sample_size, 0.75, 0.75),
    density = function(values) dbeta(values, 0.75, 0.75)
  )
)
evaluation_points <- sort(unique(c(
  10^seq(-10, -2, length.out = 101L),
  seq(0.01, 0.99, length.out = 301L),
  1 - 10^seq(-2, -10, length.out = 101L)
)))

pdf(output_file, width = 12, height = 8)
par(mfrow = c(3, 4), mar = c(4.5, 4.2, 2.4, 0.5), oma = c(0, 0, 3.5, 0))
for (sample_size in sort(unique(results$sample_size))) {
  for (distribution in names(validation_distributions)) {
    scenario_results <- results[
      results$sample_size == sample_size &
        results$distribution == distribution,
    ]
    current_result <- scenario_results[scenario_results$method == "current", ]
    selected_result <- scenario_results[scenario_results$method == "selected", ]
    ratios <- c(
      ISE = selected_result$global_ise / current_result$global_ise,
      bias = selected_result$global_squared_bias /
        current_result$global_squared_bias,
      variance = selected_result$global_variance /
        current_result$global_variance
    )
    barplot(
      ratios, col = c("grey45", "#E69F00", "#56B4E9"),
      ylim = range(c(0, 1, ratios)), ylab = "selected / current",
      main = sprintf("%s, n = %d", distribution, sample_size)
    )
    abline(h = 1, lty = 2)
  }
}
mtext(
  sprintf("Selected estimator: %d paired replications", replication_count),
  outer = TRUE, font = 2, cex = 1.1
)

par(mfrow = c(3, 4), mar = c(4.5, 4.2, 2.4, 0.5), oma = c(0, 0, 3.5, 0))
for (sample_size in sort(unique(results$sample_size))) {
  for (distribution in names(validation_distributions)) {
    selected_result <- results[
      results$sample_size == sample_size &
        results$distribution == distribution &
        results$method == "selected",
    ]
    barplot(
      rbind(
        bulk = c(selected_result$lower_bulk_rate, selected_result$upper_bulk_rate),
        finite = c(
          selected_result$lower_finite_rate,
          selected_result$upper_finite_rate
        )
      ),
      names.arg = c("lower", "upper"), col = c("grey60", "#0072B2"),
      ylim = c(0, 1), ylab = "selection rate",
      main = sprintf("%s, n = %d", distribution, sample_size)
    )
    if (sample_size == min(results$sample_size) &&
        distribution == names(validation_distributions)[[1L]]) {
      legend(
        "topright", c("bulk", "finite"),
        fill = c("grey60", "#0072B2"), bty = "n", cex = 0.75
      )
    }
  }
}
mtext("Effective endpoint components", outer = TRUE, font = 2, cex = 1.1)

sample_size <- 1000L
for (distribution_index in seq_along(validation_distributions)) {
  true_density <- validation_distributions[[distribution_index]]$density(
    evaluation_points
  )
  par(mfrow = c(2, 3), mar = c(4, 4.2, 2.3, 0.6), oma = c(0, 0, 3, 0))
  for (replication in 1:2) {
    set.seed(
      20260906L + 100000L * distribution_index + 1000L * 2L + replication
    )
    observations <- validation_distributions[[distribution_index]]$random(
      sample_size
    )
    selected_fit <- fit_bounded_selected_estimator(
      observations,
      evaluation_points
    )
    estimates <- list(
      current = normalize_density(
        evaluation_points,
        dkde1d(
          evaluation_points,
          kde1d(observations, xmin = 0, xmax = 1)
        )
      ),
      selected = selected_fit$density
    )
    plot(
      evaluation_points, true_density, type = "l", lwd = 3, col = "black",
      ylim = c(0, quantile(c(true_density, unlist(estimates)), 0.98)),
      xlab = "x", ylab = "density",
      main = sprintf(
        "realization %d: L=%s, U=%s",
        replication,
        selected_fit$expert_choices[["lower"]],
        selected_fit$expert_choices[["upper"]]
      )
    )
    for (method in names(estimates)) {
      lines(
        evaluation_points, estimates[[method]], col = method_colors[[method]],
        lwd = if (method == "selected") 2.5 else 1.5,
        lty = match(method, names(estimates))
      )
    }
    if (replication == 1L) {
      legend(
        "topright", c("truth", "current", "selected"),
        col = c("black", method_colors), lty = c(1, 1, 2),
        lwd = c(3, 1.5, 2.5), bty = "n", cex = 0.75
      )
    }
    for (boundary in c("left", "right")) {
      boundary_points <- if (boundary == "left") {
        evaluation_points >= 1e-6 & evaluation_points <= 0.1
      } else {
        evaluation_points >= 0.9 & evaluation_points <= 1 - 1e-6
      }
      boundary_distances <- if (boundary == "left") {
        evaluation_points[boundary_points]
      } else {
        1 - evaluation_points[boundary_points]
      }
      boundary_order <- order(boundary_distances)
      plot(
        boundary_distances[boundary_order],
        true_density[boundary_points][boundary_order],
        type = "l", log = "x", lwd = 3, col = "black",
        ylim = c(0, quantile(c(
          true_density[boundary_points],
          unlist(lapply(estimates, function(density) density[boundary_points]))
        ), 0.99)),
        xlab = "distance from boundary", ylab = "density", main = boundary
      )
      for (method in names(estimates)) {
        lines(
          boundary_distances[boundary_order],
          estimates[[method]][boundary_points][boundary_order],
          col = method_colors[[method]],
          lwd = if (method == "selected") 2.5 else 1.5,
          lty = match(method, names(estimates))
        )
      }
    }
  }
  mtext(
    sprintf(
      "%s, n = %d; effective choices shown in panel titles",
      names(validation_distributions)[[distribution_index]], sample_size
    ),
    outer = TRUE, font = 2, cex = 1.1
  )
}
dev.off()
message("wrote ", normalizePath(output_file))
