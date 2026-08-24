# Plot the paired endpoint-classifier ablation and representative fits.
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
  "endpoint-classifier-ablation.csv"
)
output_file <- if (length(arguments) >= 2L) arguments[[2L]] else paste0(
  "inst/include/kde1d-cpp/diagnostic/boundary-repair/plots/",
  "endpoint-classifier-ablation.pdf"
)
if (length(arguments) > 2L) {
  stop("usage: plot_endpoint_classifier_ablation.R [input.csv] [output.pdf]")
}

results <- read.csv(input_file)
classifier_methods <- c(
  "single-95", "single-1.5-95", "single-2-95",
  "median-90", "vote-90", "vote-95"
)
method_labels <- c(
  `single-95` = "1x", `single-1.5-95` = "1.5x", `single-2-95` = "2x",
  `median-90` = "median", `vote-90` = "vote 90", `vote-95` = "vote 95"
)
method_colors <- c(
  `single-95` = "#009E73", `single-1.5-95` = "#56B4E9",
  `single-2-95` = "#E69F00", `median-90` = "#CC79A7",
  `vote-90` = "#7B3294", `vote-95` = "#D55E00"
)
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

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
pdf(output_file, width = 12, height = 8)

par(mfrow = c(3, 4), mar = c(4.5, 4.2, 2.4, 0.5), oma = c(0, 0, 3.5, 0))
for (sample_size in sort(unique(results$sample_size))) {
  for (distribution in names(validation_distributions)) {
    scenario_results <- results[
      results$sample_size == sample_size &
        results$distribution == distribution,
    ]
    current_ise <- scenario_results$global_ise[
      scenario_results$method == "kde-current"
    ]
    scenario_results <- scenario_results[
      match(classifier_methods, scenario_results$method),
    ]
    ratios <- scenario_results$global_ise / current_ise
    barplot(
      ratios,
      names.arg = method_labels[classifier_methods],
      col = method_colors[classifier_methods], las = 2, cex.names = 0.7,
      ylim = range(c(0, 1, ratios)), ylab = "global ISE / current",
      main = sprintf("%s, n = %d", distribution, sample_size)
    )
    abline(h = 1, lty = 2)
  }
}
mtext(
  sprintf(
    "Classifier ISE: %d paired replications",
    unique(results$replications)
  ),
  outer = TRUE, font = 2, cex = 1.1
)

par(mfrow = c(3, 4), mar = c(4.5, 4.2, 2.4, 0.5), oma = c(0, 0, 3.5, 0))
for (sample_size in sort(unique(results$sample_size))) {
  for (distribution in names(validation_distributions)) {
    scenario_results <- results[
      results$sample_size == sample_size &
        results$distribution == distribution &
        results$method %in% classifier_methods,
    ]
    scenario_results <- scenario_results[
      match(classifier_methods, scenario_results$method),
    ]
    barplot(
      rowMeans(scenario_results[, c("lower_accuracy", "upper_accuracy")]),
      names.arg = method_labels[classifier_methods],
      col = method_colors[classifier_methods], las = 2, cex.names = 0.7,
      ylim = c(0, 1), ylab = "mean endpoint accuracy",
      main = sprintf("%s, n = %d", distribution, sample_size)
    )
  }
}
mtext(
  "Classification accuracy (finite versus bulk endpoint component)",
  outer = TRUE, font = 2, cex = 1.1
)

visual_methods <- c("single-95", "single-2-95", "vote-90")
sample_size <- 1000L
for (distribution_index in seq_along(validation_distributions)) {
  true_density <- validation_distributions[[distribution_index]]$density(
    evaluation_points
  )
  par(mfrow = c(2, 3), mar = c(4, 4.2, 2.4, 0.6), oma = c(0, 0, 3, 0))
  for (replication in 1:2) {
    set.seed(
      20260909L + 100000L * distribution_index + 1000L * 2L + replication
    )
    observations <- validation_distributions[[distribution_index]]$random(
      sample_size
    )
    estimator_components <- fit_bounded_estimator_components(
      observations,
      evaluation_points
    )
    selections <- list(
      `single-95` = select_bounded_endpoint_experts_single_scale(observations),
      `single-2-95` = select_bounded_endpoint_experts_single_scale(
        observations,
        tail_multiplier = 2
      ),
      `vote-90` = select_bounded_endpoint_experts_multiscale(
        observations,
        evidence_threshold = qnorm(0.90),
        rule = "vote"
      )
    )
    estimates <- list(
      current = normalize_density(
        evaluation_points,
        dkde1d(
          evaluation_points,
          kde1d(observations, xmin = 0, xmax = 1)
        )
      )
    )
    for (method in visual_methods) {
      estimator_components$selection <- selections[[method]]
      estimates[[method]] <- fuse_bounded_estimator_components(
        evaluation_points,
        estimator_components
      )$density
    }
    plot(
      evaluation_points, true_density, type = "l", lwd = 3, col = "black",
      ylim = c(0, quantile(c(true_density, unlist(estimates)), 0.98)),
      xlab = "x", ylab = "density",
      main = sprintf("realization %d", replication)
    )
    lines(evaluation_points, estimates$current, col = "grey45", lwd = 1.5)
    for (method in visual_methods) {
      lines(
        evaluation_points, estimates[[method]], col = method_colors[[method]],
        lwd = if (method == "single-2-95") 2.5 else 1.5,
        lty = match(method, visual_methods)
      )
    }
    if (replication == 1L) {
      legend(
        "topright", c("truth", "current", method_labels[visual_methods]),
        col = c("black", "grey45", method_colors[visual_methods]),
        lty = c(1, 1, seq_along(visual_methods)),
        lwd = c(3, 1.5, 1.5, 2.5, 1.5), bty = "n", cex = 0.72
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
        ), 0.98)),
        xlab = "distance from boundary", ylab = "density", main = boundary
      )
      lines(
        boundary_distances[boundary_order],
        estimates$current[boundary_points][boundary_order],
        col = "grey45", lwd = 1.5
      )
      for (method in visual_methods) {
        lines(
          boundary_distances[boundary_order],
          estimates[[method]][boundary_points][boundary_order],
          col = method_colors[[method]],
          lwd = if (method == "single-2-95") 2.5 else 1.5,
          lty = match(method, visual_methods)
        )
      }
    }
  }
  mtext(
    sprintf("%s, n = %d", names(validation_distributions)[distribution_index], sample_size),
    outer = TRUE, font = 2, cex = 1.1
  )
}

dev.off()
message("wrote ", normalizePath(output_file))
