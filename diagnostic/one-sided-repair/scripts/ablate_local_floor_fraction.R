# Paired ablation of the local probability mass used for the bandwidth floor.
# Run from the root of the kde1d R package.

if (!requireNamespace("pkgload", quietly = TRUE) || !file.exists("DESCRIPTION")) {
  stop("run from the kde1d package root with pkgload installed")
}
pkgload::load_all(".", quiet = TRUE)
source(paste0(
  "inst/include/kde1d-cpp/diagnostic/one-sided-repair/scripts/",
  "one_sided_estimator.R"
))

arguments <- commandArgs(trailingOnly = TRUE)
replication_count <- if (length(arguments) >= 1L) {
  as.integer(arguments[[1L]])
} else {
  30L
}
output_file <- if (length(arguments) >= 2L) arguments[[2L]] else paste0(
  "inst/include/kde1d-cpp/diagnostic/one-sided-repair/results/",
  "local-floor-fraction-ablation.csv"
)
if (length(arguments) > 2L || !is.finite(replication_count) ||
    replication_count < 2L) {
  stop("usage: ablate_local_floor_fraction.R [replications] [output.csv]")
}

validation_distributions <- list(
  exponential = list(random = rexp, density = dexp, quantile = qexp),
  half_normal = list(
    random = function(sample_size) abs(rnorm(sample_size)),
    density = function(values) 2 * dnorm(values),
    quantile = function(probabilities) qnorm((probabilities + 1) / 2)
  ),
  gamma_2 = list(
    random = function(sample_size) rgamma(sample_size, 2),
    density = function(values) dgamma(values, 2),
    quantile = function(probabilities) qgamma(probabilities, 2)
  ),
  gamma_0.75 = list(
    random = function(sample_size) rgamma(sample_size, 0.75),
    density = function(values) dgamma(values, 0.75),
    quantile = function(probabilities) qgamma(probabilities, 0.75)
  ),
  lomax_2 = list(
    random = function(sample_size) (1 - runif(sample_size))^(-0.5) - 1,
    density = function(values) 2 / (1 + values)^3,
    quantile = function(probabilities) (1 - probabilities)^(-0.5) - 1
  )
)
sample_sizes <- c(100L, 1000L)
local_floor_fractions <- c(
  0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.75, 1
)
fraction_methods <- paste0(
  "fixed-",
  format(100 * local_floor_fractions, trim = TRUE),
  "pct"
)
candidate_methods <- c(fraction_methods, "shrinking")
methods <- c("current", candidate_methods)
evaluation_probabilities <- sort(unique(c(
  10^seq(-10, -2, length.out = 101L),
  seq(0.01, 0.999, length.out = 401L),
  0.9999, 0.99999
)))
result_rows <- list()

summarize_fraction_estimates <- function(
    evaluation_points,
    estimates,
    true_density,
    region_points) {
  mean_density <- rowMeans(estimates)
  integrated_errors <- vapply(seq_len(ncol(estimates)), function(replication) {
    one_sided_trapezoid_integral(
      evaluation_points[region_points],
      (estimates[region_points, replication] - true_density[region_points])^2
    )
  }, numeric(1L))
  c(
    ise = mean(integrated_errors),
    squared_bias = one_sided_trapezoid_integral(
      evaluation_points[region_points],
      (mean_density[region_points] - true_density[region_points])^2
    ),
    variance = one_sided_trapezoid_integral(
      evaluation_points[region_points],
      rowMeans((estimates[region_points, , drop = FALSE] -
        mean_density[region_points])^2)
    )
  )
}

for (distribution_index in seq_along(validation_distributions)) {
  evaluation_points <- validation_distributions[[distribution_index]]$quantile(
    evaluation_probabilities
  )
  true_density <- validation_distributions[[distribution_index]]$density(
    evaluation_points
  )
  boundary_points <- evaluation_probabilities <= 0.1
  for (sample_size in sample_sizes) {
    estimates <- setNames(lapply(methods, function(method) {
      matrix(NA_real_, length(evaluation_points), replication_count)
    }), methods)
    finite_choices <- logical(replication_count)
    floor_ratios <- setNames(lapply(candidate_methods, function(method) {
      rep(NA_real_, replication_count)
    }), candidate_methods)
    for (replication in seq_len(replication_count)) {
      set.seed(
        20260914L + 100000L * distribution_index +
          1000L * match(sample_size, sample_sizes) + replication
      )
      observations <- validation_distributions[[distribution_index]]$random(
        sample_size
      )
      estimates$current[, replication] <- normalize_one_sided_density(
        evaluation_points,
        dkde1d(evaluation_points, kde1d(observations, xmin = 0))
      )
      for (fraction_index in seq_along(local_floor_fractions)) {
        method <- fraction_methods[[fraction_index]]
        candidate_fit <- fit_one_sided_candidate(
          observations,
          evaluation_points,
          bandwidth_floor = "local",
          local_floor_fraction = local_floor_fractions[[fraction_index]]
        )
        estimates[[method]][, replication] <- candidate_fit$density
        floor_ratios[[method]][[replication]] <- mean(
          candidate_fit$boundary_bandwidths /
            candidate_fit$selected_boundary_bandwidths
        )
        if (fraction_index == 1L) {
          finite_choices[[replication]] <-
            candidate_fit$expert_choice == "finite"
        }
      }
      shrinking_fit <- fit_one_sided_candidate(
        observations,
        evaluation_points,
        bandwidth_floor = "shrinking"
      )
      estimates$shrinking[, replication] <- shrinking_fit$density
      floor_ratios$shrinking[[replication]] <- mean(
        shrinking_fit$boundary_bandwidths /
          shrinking_fit$selected_boundary_bandwidths
      )
    }
    for (method in methods) {
      global_summary <- summarize_fraction_estimates(
        evaluation_points,
        estimates[[method]],
        true_density,
        rep(TRUE, length(evaluation_points))
      )
      boundary_summary <- summarize_fraction_estimates(
        evaluation_points,
        estimates[[method]],
        true_density,
        boundary_points
      )
      result_rows[[length(result_rows) + 1L]] <- data.frame(
        distribution = names(validation_distributions)[[distribution_index]],
        sample_size = sample_size,
        replications = replication_count,
        method = method,
        global_ise = global_summary[["ise"]],
        global_squared_bias = global_summary[["squared_bias"]],
        global_variance = global_summary[["variance"]],
        boundary_ise = boundary_summary[["ise"]],
        boundary_squared_bias = boundary_summary[["squared_bias"]],
        boundary_variance = boundary_summary[["variance"]],
        finite_rate = if (method %in% candidate_methods) {
          mean(finite_choices)
        } else {
          NA_real_
        },
        mean_floor_ratio = if (
          method %in% candidate_methods &&
          any(is.finite(floor_ratios[[method]]))
        ) {
          mean(floor_ratios[[method]], na.rm = TRUE)
        } else {
          NA_real_
        },
        minimum_density = min(estimates[[method]]),
        maximum_mass_error = max(abs(vapply(
          seq_len(replication_count),
          function(replication) {
            one_sided_trapezoid_integral(
              evaluation_points,
              estimates[[method]][, replication]
            ) - 1
          },
          numeric(1L)
        )))
      )
    }
    message(
      names(validation_distributions)[[distribution_index]],
      ", n = ", sample_size
    )
  }
}

results <- do.call(rbind, result_rows)
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write.csv(results, output_file, row.names = FALSE)
message("wrote ", normalizePath(output_file))
