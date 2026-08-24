# Evaluate one native kde1d revision. Run once at the baseline commit and once
# at the candidate commit, then compare rows by scenario, size, and replication.

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
  points <- sort(unique(c(
    1e-10,
    scenario$quantile(10^seq(-8, -1, length.out = 100L)),
    scenario$quantile(seq(0.1, 0.9999, length.out = 401L))
  )))
  truth <- scenario$density(points)
  boundary <- points <= scenario$quantile(0.1)
  for (sample_size in sample_sizes) {
    estimates <- matrix(NA_real_, length(points), replications)
    for (replication in seq_len(replications)) {
      set.seed(20260824L + 100000L * scenario_index +
                 1000L * match(sample_size, sample_sizes) + replication)
      fit <- kde1d(scenario$random(sample_size), xmin = 0)
      estimates[, replication] <- dkde1d(points, fit)
      result_rows[[length(result_rows) + 1L]] <- data.frame(
        support = "one-sided",
        scenario = names(one_sided)[[scenario_index]],
        sample_size = sample_size,
        replication = replication,
        method = method,
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
      scenario = names(one_sided)[[scenario_index]],
      sample_size = sample_size,
      method = method,
      as.list(summarize_estimates(points, estimates, truth, boundary))
    )
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
  for (sample_size in sample_sizes) {
    estimates <- matrix(NA_real_, length(points), replications)
    for (replication in seq_len(replications)) {
      set.seed(20260825L + 100000L * scenario_index +
                 1000L * match(sample_size, sample_sizes) + replication)
      fit <- kde1d(
        rbeta(sample_size, shapes[[1L]], shapes[[2L]]),
        xmin = 0,
        xmax = 1
      )
      estimates[, replication] <- dkde1d(points, fit)
      result_rows[[length(result_rows) + 1L]] <- data.frame(
        support = "two-sided",
        scenario = names(two_sided)[[scenario_index]],
        sample_size = sample_size,
        replication = replication,
        method = method,
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
      scenario = names(two_sided)[[scenario_index]],
      sample_size = sample_size,
      method = method,
      as.list(summarize_estimates(points, estimates, truth, boundary))
    )
  }
}

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write.csv(do.call(rbind, result_rows), output_file, row.names = FALSE)
write.csv(
  do.call(rbind, summary_rows),
  sub("\\.csv$", "-summary.csv", output_file),
  row.names = FALSE
)
