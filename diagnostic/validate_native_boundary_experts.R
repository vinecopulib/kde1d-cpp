# Evaluate one native kde1d method. Install the pre-PR revision for `pre-pr`;
# install the candidate revision for `bulk` and `expert`. The pre-PR one-sided
# fit is the production log-transform estimator, not a local emulation of it.

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) < 2L || length(arguments) > 4L) {
  stop("usage: script method output.csv [replications] [library-path]")
}
if (length(arguments) == 4L) {
  .libPaths(c(arguments[[4L]], .libPaths()))
}
library(kde1d)

method <- arguments[[1L]]
output_file <- arguments[[2L]]
if (!method %in% c("pre-pr", "bulk", "expert")) {
  stop("method must be one of: pre-pr, bulk, expert")
}
replications <- if (length(arguments) >= 3L) {
  as.integer(arguments[[3L]])
} else {
  100L
}
if (!is.finite(replications) || replications < 2L) {
  stop("replications must be at least two")
}

integrate_values <- function(points, values) {
  sum(diff(points) * (values[-length(values)] + values[-1L]) / 2)
}

fit_kde1d <- function(observations, degree,
                      xmin = NA_real_, xmax = NA_real_) {
  arguments <- list(x = observations, xmin = xmin, xmax = xmax, deg = degree)
  if (method != "pre-pr") {
    arguments$boundary_repair <- method == "expert"
  }
  do.call(kde1d, arguments)
}

classify_finite_endpoint <- function(distances) {
  distances <- sort(distances)
  k <- min(length(distances) - 1L, ceiling(2 * sqrt(length(distances))))
  reference_distance <- distances[[k + 1L]]
  if (!is.finite(reference_distance) || reference_distance <= 0) {
    return(FALSE)
  }
  denominator <- sum(log(
    reference_distance /
      pmax(distances[seq_len(k)], .Machine$double.eps * reference_distance)
  ))
  if (!is.finite(denominator) || denominator <= 0) {
    return(FALSE)
  }
  beta <- k / denominator
  beta >= 0.9 && beta * (1 - qnorm(0.95) / sqrt(k)) <= 1
}

one_sided <- list(
  exponential = list(
    random = function(n) rexp(n),
    density = dexp,
    quantile = qexp
  ),
  half_normal = list(
    random = function(n) abs(rnorm(n)),
    density = function(x) 2 * dnorm(x),
    quantile = function(p) qnorm((p + 1) / 2)
  ),
  gamma_2 = list(
    random = function(n) rgamma(n, 2),
    density = function(x) dgamma(x, 2),
    quantile = function(p) qgamma(p, 2)
  ),
  gamma_075 = list(
    random = function(n) rgamma(n, 0.75),
    density = function(x) dgamma(x, 0.75),
    quantile = function(p) qgamma(p, 0.75)
  ),
  lomax_2 = list(
    random = function(n) (1 - runif(n))^(-0.5) - 1,
    density = function(x) 2 / (1 + x)^3,
    quantile = function(p) (1 - p)^(-0.5) - 1
  )
)
two_sided <- list(
  uniform = c(1, 1),
  beta_1_2 = c(1, 2),
  beta_2_1 = c(2, 1),
  beta_075_2 = c(0.75, 2),
  beta_2_075 = c(2, 0.75),
  beta_075_075 = c(0.75, 0.75),
  beta_2_2 = c(2, 2)
)
sample_sizes <- c(25L, 100L, 1000L, 2000L)
degrees <- 0:2
one_sided_scales <- c(1e-4, 1, 1e4)
one_sided_directions <- c("lower", "upper")
result_rows <- list()
summary_rows <- list()

summarize_estimates <- function(points, estimates, truth, boundary) {
  mean_density <- rowMeans(estimates)
  c(
    global_squared_bias = integrate_values(
      points, (mean_density - truth)^2
    ),
    global_variance = integrate_values(
      points, rowMeans((estimates - mean_density)^2)
    ),
    boundary_squared_bias = integrate_values(
      points[boundary], (mean_density[boundary] - truth[boundary])^2
    ),
    boundary_variance = integrate_values(
      points[boundary],
      rowMeans((estimates[boundary, , drop = FALSE] -
                  mean_density[boundary])^2)
    )
  )
}

for (scenario_index in seq_along(one_sided)) {
  scenario <- one_sided[[scenario_index]]
  base_points <- sort(unique(c(
    1e-10,
    scenario$quantile(10^seq(-8, -1, length.out = 100L)),
    scenario$quantile(seq(0.1, 0.9999, length.out = 401L))
  )))
  for (scale_index in seq_along(one_sided_scales)) {
    scale <- one_sided_scales[[scale_index]]
    scaled_points <- scale * base_points
    scaled_truth <- scenario$density(base_points) / scale
    for (direction_index in seq_along(one_sided_directions)) {
      direction <- one_sided_directions[[direction_index]]
      if (direction == "lower") {
        points <- scaled_points
        truth <- scaled_truth
        boundary <- points <= scale * scenario$quantile(0.1)
      } else {
        points <- -rev(scaled_points)
        truth <- rev(scaled_truth)
        boundary <- points >= -scale * scenario$quantile(0.1)
      }
      for (degree in degrees) {
        for (sample_size in sample_sizes) {
          estimates <- matrix(NA_real_, length(points), replications)
          for (replication in seq_len(replications)) {
            # Reuse the base sample across scales and reflections so departures
            # from scale/reflection equivariance are paired directly.
            set.seed(20260824L + 10000000L * scenario_index +
                       1000L * match(sample_size, sample_sizes) + replication)
            observations <- scale * scenario$random(sample_size)
            repair_endpoint <- classify_finite_endpoint(observations)
            if (direction == "lower") {
              fit <- fit_kde1d(observations, degree, xmin = 0)
            } else {
              fit <- fit_kde1d(-observations, degree, xmax = 0)
            }
            estimates[, replication] <- dkde1d(points, fit)
            result_rows[[length(result_rows) + 1L]] <- data.frame(
              support = "one-sided",
              direction = direction,
              scale = scale,
              scenario = names(one_sided)[[scenario_index]],
              degree = degree,
              sample_size = sample_size,
              replication = replication,
              method = method,
              repair_lower = direction == "lower" && repair_endpoint,
              repair_upper = direction == "upper" && repair_endpoint,
              global_ise = integrate_values(
                points, (estimates[, replication] - truth)^2
              ),
              boundary_ise = integrate_values(
                points[boundary],
                (estimates[boundary, replication] - truth[boundary])^2
              ),
              edf = fit$edf,
              loglik = fit$loglik
            )
          }
          summary_rows[[length(summary_rows) + 1L]] <- data.frame(
            support = "one-sided",
            direction = direction,
            scale = scale,
            scenario = names(one_sided)[[scenario_index]],
            degree = degree,
            sample_size = sample_size,
            method = method,
            as.list(summarize_estimates(points, estimates, truth, boundary))
          )
        }
      }
    }
  }
}

points <- sort(unique(c(
  10^seq(-10, -2, length.out = 101L),
  seq(0.01, 0.99, length.out = 301L),
  1 - 10^seq(-2, -10, length.out = 101L)
)))
boundary <- points <= 0.1 | points >= 0.9
for (scenario_index in seq_along(two_sided)) {
  shapes <- two_sided[[scenario_index]]
  truth <- dbeta(points, shapes[[1L]], shapes[[2L]])
  for (degree in degrees) {
    for (sample_size in sample_sizes) {
      estimates <- matrix(NA_real_, length(points), replications)
      for (replication in seq_len(replications)) {
        set.seed(20260825L + 100000L * scenario_index +
                   1000L * match(sample_size, sample_sizes) + replication)
        observations <- rbeta(sample_size, shapes[[1L]], shapes[[2L]])
        fit <- fit_kde1d(
          observations,
          degree,
          xmin = 0,
          xmax = 1
        )
        estimates[, replication] <- dkde1d(points, fit)
        result_rows[[length(result_rows) + 1L]] <- data.frame(
          support = "two-sided",
          direction = NA_character_,
          scale = 1,
          scenario = names(two_sided)[[scenario_index]],
          degree = degree,
          sample_size = sample_size,
          replication = replication,
          method = method,
          repair_lower = classify_finite_endpoint(observations),
          repair_upper = classify_finite_endpoint(1 - observations),
          global_ise = integrate_values(
            points, (estimates[, replication] - truth)^2
          ),
          boundary_ise = integrate_values(
            points[boundary],
            (estimates[boundary, replication] - truth[boundary])^2
          ),
          edf = fit$edf,
          loglik = fit$loglik
        )
      }
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        support = "two-sided",
        direction = NA_character_,
        scale = 1,
        scenario = names(two_sided)[[scenario_index]],
        degree = degree,
        sample_size = sample_size,
        method = method,
        as.list(summarize_estimates(points, estimates, truth, boundary))
      )
    }
  }
}

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write.csv(do.call(rbind, result_rows), output_file, row.names = FALSE)
write.csv(
  do.call(rbind, summary_rows),
  sub("\\.csv$", "-summary.csv", output_file),
  row.names = FALSE
)
