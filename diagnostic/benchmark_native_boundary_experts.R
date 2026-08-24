arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 3L) {
  stop("usage: script method output.csv library-path")
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

make_observations <- function(case, n) {
  probabilities <- seq(0.5 / n, 1 - 0.5 / n, length.out = n)
  switch(
    case,
    one_active = -log1p(-probabilities),
    one_fallback = sqrt(probabilities),
    two_active = probabilities,
    two_fallback = qbeta(probabilities, 2, 2)
  )
}

rows <- list()
for (case in c("one_active", "one_fallback", "two_active", "two_fallback")) {
  for (n in c(100L, 1000L, 10000L)) {
    observations <- make_observations(case, n)
    one_sided <- startsWith(case, "one_")
    xmin <- 0
    xmax <- if (one_sided) NaN else 1
    iterations <- if (n == 100L) 500L else if (n == 1000L) 200L else 30L
    fit_kde1d(observations, xmin, xmax)
    for (batch in seq_len(7L)) {
      elapsed <- system.time({
        for (iteration in seq_len(iterations)) {
          fit_kde1d(observations, xmin, xmax)
        }
      })[["elapsed"]]
      rows[[length(rows) + 1L]] <- data.frame(
        method = method,
        case = case,
        sample_size = n,
        batch = batch,
        milliseconds = 1000 * elapsed / iterations
      )
    }
  }
}

dir.create(dirname(arguments[[2L]]), recursive = TRUE, showWarnings = FALSE)
write.csv(do.call(rbind, rows), arguments[[2L]], row.names = FALSE)
