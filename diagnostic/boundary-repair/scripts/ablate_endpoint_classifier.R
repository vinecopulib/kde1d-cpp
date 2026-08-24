# Paired ablation of order-statistic endpoint classifiers.
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
replication_count <- if (length(arguments) >= 1L) {
  as.integer(arguments[[1L]])
} else {
  100L
}
output_file <- if (length(arguments) >= 2L) arguments[[2L]] else paste0(
  "inst/include/kde1d-cpp/diagnostic/boundary-repair/results/",
  "endpoint-classifier-ablation.csv"
)
if (length(arguments) > 2L || !is.finite(replication_count) ||
    replication_count < 2L) {
  stop("usage: ablate_endpoint_classifier.R [replications] [output.csv]")
}

make_beta_distribution <- function(shape1, shape2, endpoint_truth) {
  list(
    random = function(sample_size) rbeta(sample_size, shape1, shape2),
    density = function(values) dbeta(values, shape1, shape2),
    endpoint_truth = endpoint_truth
  )
}
validation_distributions <- list(
  uniform = list(
    random = runif,
    density = dunif,
    endpoint_truth = c(lower = "finite", upper = "finite")
  ),
  `beta-1-2` = make_beta_distribution(
    1, 2, c(lower = "finite", upper = "bulk")
  ),
  `beta-0.75-2` = make_beta_distribution(
    0.75, 2, c(lower = "bulk", upper = "bulk")
  ),
  `beta-0.75-0.75` = make_beta_distribution(
    0.75, 0.75, c(lower = "bulk", upper = "bulk")
  )
)
sample_sizes <- c(100L, 1000L, 2000L)
classifier_methods <- c(
  "single-95", "single-1.5-95", "single-2-95",
  "median-90", "vote-90", "vote-95"
)
methods <- c("kde-current", classifier_methods)
evaluation_points <- sort(unique(c(
  10^seq(-10, -2, length.out = 101L),
  seq(0.01, 0.99, length.out = 301L),
  1 - 10^seq(-2, -10, length.out = 101L)
)))
left_points <- evaluation_points <= 0.1
right_points <- evaluation_points >= 0.9
result_rows <- list()

make_classifier_selections <- function(observations) {
  list(
    `single-95` = select_bounded_endpoint_experts_single_scale(observations),
    `single-1.5-95` = select_bounded_endpoint_experts_single_scale(
      observations,
      tail_multiplier = 1.5
    ),
    `single-2-95` = select_bounded_endpoint_experts_single_scale(
      observations,
      tail_multiplier = 2
    ),
    `median-90` = select_bounded_endpoint_experts_multiscale(
      observations,
      evidence_threshold = qnorm(0.90),
      rule = "median"
    ),
    `vote-90` = select_bounded_endpoint_experts_multiscale(
      observations,
      evidence_threshold = qnorm(0.90),
      rule = "vote"
    ),
    `vote-95` = select_bounded_endpoint_experts_multiscale(
      observations,
      evidence_threshold = qnorm(0.95),
      rule = "vote"
    )
  )
}

summarize_estimates <- function(estimates, true_density, region_points) {
  mean_density <- rowMeans(estimates)
  integrated_errors <- vapply(
    seq_len(ncol(estimates)),
    function(replication) {
      trapezoid_integral(
        evaluation_points[region_points],
        (estimates[region_points, replication] -
          true_density[region_points])^2
      )
    },
    numeric(1L)
  )
  c(
    ise = mean(integrated_errors),
    squared_bias = trapezoid_integral(
      evaluation_points[region_points],
      (mean_density[region_points] - true_density[region_points])^2
    ),
    variance = trapezoid_integral(
      evaluation_points[region_points],
      rowMeans((estimates[region_points, , drop = FALSE] -
        mean_density[region_points])^2)
    )
  )
}

for (distribution_index in seq_along(validation_distributions)) {
  true_density <- validation_distributions[[distribution_index]]$density(
    evaluation_points
  )
  for (sample_size in sample_sizes) {
    estimates <- setNames(
      lapply(methods, function(method) {
        matrix(NA_real_, length(evaluation_points), replication_count)
      }),
      methods
    )
    choices <- setNames(
      lapply(classifier_methods, function(method) {
        matrix(
          NA_character_, replication_count, 2L,
          dimnames = list(NULL, c("lower", "upper"))
        )
      }),
      classifier_methods
    )
    for (replication in seq_len(replication_count)) {
      set.seed(
        20260909L + 100000L * distribution_index +
          1000L * match(sample_size, sample_sizes) + replication
      )
      observations <- validation_distributions[[distribution_index]]$random(
        sample_size
      )
      estimates$`kde-current`[, replication] <- normalize_density(
        evaluation_points,
        dkde1d(
          evaluation_points,
          kde1d(observations, xmin = 0, xmax = 1)
        )
      )
      estimator_components <- fit_bounded_estimator_components(
        observations,
        evaluation_points
      )
      classifier_selections <- make_classifier_selections(observations)
      for (method in classifier_methods) {
        estimator_components$selection <- classifier_selections[[method]]
        estimates[[method]][, replication] <-
          fuse_bounded_estimator_components(
            evaluation_points,
            estimator_components
          )$density
        choices[[method]][replication, ] <-
          classifier_selections[[method]]$expert_choices
      }
    }
    for (method in methods) {
      global_summary <- summarize_estimates(
        estimates[[method]], true_density, rep(TRUE, length(evaluation_points))
      )
      left_summary <- summarize_estimates(
        estimates[[method]], true_density, left_points
      )
      right_summary <- summarize_estimates(
        estimates[[method]], true_density, right_points
      )
      endpoint_truth <- validation_distributions[[distribution_index]]$endpoint_truth
      result_rows[[length(result_rows) + 1L]] <- data.frame(
        distribution = names(validation_distributions)[[distribution_index]],
        sample_size = sample_size,
        replications = replication_count,
        method = method,
        global_ise = global_summary[["ise"]],
        global_squared_bias = global_summary[["squared_bias"]],
        global_variance = global_summary[["variance"]],
        left_ise = left_summary[["ise"]],
        left_squared_bias = left_summary[["squared_bias"]],
        left_variance = left_summary[["variance"]],
        right_ise = right_summary[["ise"]],
        right_squared_bias = right_summary[["squared_bias"]],
        right_variance = right_summary[["variance"]],
        lower_accuracy = if (method == "kde-current") {
          NA_real_
        } else {
          mean(choices[[method]][, "lower"] == endpoint_truth[["lower"]])
        },
        upper_accuracy = if (method == "kde-current") {
          NA_real_
        } else {
          mean(choices[[method]][, "upper"] == endpoint_truth[["upper"]])
        },
        lower_finite_rate = if (method == "kde-current") NA_real_ else
          mean(choices[[method]][, "lower"] == "finite"),
        upper_finite_rate = if (method == "kde-current") NA_real_ else
          mean(choices[[method]][, "upper"] == "finite"),
        minimum_density = min(estimates[[method]]),
        maximum_mass_error = max(abs(vapply(
          seq_len(replication_count),
          function(replication) {
            trapezoid_integral(
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
