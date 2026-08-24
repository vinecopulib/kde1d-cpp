arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 4L) {
  stop("usage: script method output.csv library-path selections.csv")
}
.libPaths(c(arguments[[3L]], .libPaths()))
library(kde1d)

method <- arguments[[1L]]
if (!method %in% c("pre-pr", "bulk", "expert")) {
  stop("method must be one of: pre-pr, bulk, expert")
}

fit_kde1d <- function(observations, xmin, xmax) {
  fit_arguments <- list(x = observations, xmin = xmin, xmax = xmax, deg = 2)
  if (method != "pre-pr") {
    fit_arguments$boundary_repair <- method == "expert"
  }
  do.call(kde1d, fit_arguments)
}

one_sided <- list(
  exponential = list(random = rexp, density = dexp, quantile = qexp),
  half_normal = list(
    random = function(n) abs(rnorm(n)),
    density = function(x) 2 * dnorm(x),
    quantile = function(p) qnorm((p + 1) / 2)
  ),
  gamma_075 = list(
    random = function(n) rgamma(n, 0.75),
    density = function(x) dgamma(x, 0.75),
    quantile = function(p) qgamma(p, 0.75)
  ),
  gamma_2 = list(
    random = function(n) rgamma(n, 2),
    density = function(x) dgamma(x, 2),
    quantile = function(p) qgamma(p, 2)
  )
)
two_sided <- list(
  uniform = list(shape1 = 1, shape2 = 1),
  beta_1_2 = list(shape1 = 1, shape2 = 2),
  beta_075_2 = list(shape1 = 0.75, shape2 = 2),
  beta_2_2 = list(shape1 = 2, shape2 = 2)
)
selections <- read.csv(arguments[[4L]], stringsAsFactors = FALSE)

rows <- list()
for (scenario_index in seq_along(one_sided)) {
  scenario <- one_sided[[scenario_index]]
  selection <- selections[
    selections$support == "one-sided" &
      selections$scenario == names(one_sided)[[scenario_index]],
  ]
  sample_size_index <- match(selection$sample_size, c(25L, 100L, 1000L, 2000L))
  set.seed(20260824L + 10000000L * scenario_index +
             1000L * sample_size_index + selection$replication)
  observations <- scenario$random(selection$sample_size)
  points <- sort(unique(c(
    1e-10,
    scenario$quantile(10^seq(-8, -1, length.out = 120L)),
    scenario$quantile(seq(0.1, 0.999, length.out = 300L))
  )))
  if (is.finite(scenario$density(0))) {
    points <- c(0, points)
  }
  fit <- fit_kde1d(observations, 0, NaN)
  rows[[length(rows) + 1L]] <- data.frame(
    method = method,
    support = "one-sided",
    scenario = paste0(
      names(one_sided)[[scenario_index]], "_n", selection$sample_size
    ),
    sample_size = selection$sample_size,
    x = points,
    truth = scenario$density(points),
    estimate = dkde1d(points, fit)
  )
}

points <- sort(unique(c(
  10^seq(-10, -2, length.out = 121L),
  seq(0.01, 0.99, length.out = 301L),
  1 - 10^seq(-2, -10, length.out = 121L)
)))
for (scenario_index in seq_along(two_sided)) {
  scenario <- two_sided[[scenario_index]]
  selection <- selections[
    selections$support == "two-sided" &
      selections$scenario == names(two_sided)[[scenario_index]],
  ]
  sample_size_index <- match(selection$sample_size, c(25L, 100L, 1000L, 2000L))
  set.seed(20260825L + 100000L * scenario_index +
             1000L * sample_size_index + selection$replication)
  observations <- rbeta(
    selection$sample_size, scenario$shape1, scenario$shape2
  )
  scenario_points <- points
  if (scenario$shape1 >= 1 && scenario$shape2 >= 1) {
    scenario_points <- c(0, points, 1)
  }
  fit <- fit_kde1d(observations, 0, 1)
  rows[[length(rows) + 1L]] <- data.frame(
    method = method,
    support = "two-sided",
    scenario = paste0(
      names(two_sided)[[scenario_index]], "_n", selection$sample_size
    ),
    sample_size = selection$sample_size,
    x = scenario_points,
    truth = dbeta(scenario_points, scenario$shape1, scenario$shape2),
    estimate = dkde1d(scenario_points, fit)
  )
}

dir.create(dirname(arguments[[2L]]), recursive = TRUE, showWarnings = FALSE)
write.csv(do.call(rbind, rows), arguments[[2L]], row.names = FALSE)
