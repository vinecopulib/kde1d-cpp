arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 4L) {
  stop("usage: script pre-pr.csv bulk.csv expert.csv output-prefix")
}

read_results <- function(path) {
  results <- read.csv(path, stringsAsFactors = FALSE)
  results$direction[is.na(results$direction)] <- "none"
  results
}

pre_pr <- read_results(arguments[[1L]])
bulk <- read_results(arguments[[2L]])
expert <- read_results(arguments[[3L]])
output_prefix <- arguments[[4L]]
keys <- c(
  "support", "direction", "scale", "scenario", "degree", "sample_size",
  "replication"
)

paired_ratios <- function(reference, comparison) {
  paired <- merge(
    expert,
    reference,
    by = keys,
    suffixes = c("_expert", "_reference")
  )
  groups <- split(
    seq_len(nrow(paired)),
    interaction(paired[c(
      "support", "direction", "scale", "scenario", "degree", "sample_size"
    )], drop = TRUE)
  )
  rows <- lapply(groups, function(index) {
    lapply(c("global_ise", "boundary_ise"), function(metric) {
      candidate <- paired[[paste0(metric, "_expert")]][index]
      baseline <- paired[[paste0(metric, "_reference")]][index]
      ratio <- mean(candidate) / mean(baseline)
      data.frame(
        comparison = comparison,
        metric = metric,
        paired[index[[1L]], c(
          "support", "direction", "scale", "scenario", "degree",
          "sample_size"
        )],
        ratio = ratio,
        mcse = sd(candidate - ratio * baseline) /
          (sqrt(length(index)) * mean(baseline))
      )
    })
  })
  do.call(rbind, unlist(rows, recursive = FALSE))
}

cell_ratios <- rbind(
  paired_ratios(pre_pr, "pre-pr"),
  paired_ratios(bulk, "bulk")
)

# Average scale/reflection cells within a family while retaining their paired
# covariance. Different families use independent simulation streams.
family_paired_ratios <- function(reference, comparison) {
  paired <- merge(
    expert,
    reference,
    by = keys,
    suffixes = c("_expert", "_reference")
  )
  cell_names <- c(
    "support", "direction", "scale", "scenario", "degree", "sample_size"
  )
  family_names <- c("support", "scenario", "degree", "sample_size")
  rows <- lapply(c("global_ise", "boundary_ise"), function(metric) {
    paired$cell <- interaction(paired[cell_names], drop = TRUE)
    paired$influence <- NA_real_
    paired$cell_ratio <- NA_real_
    for (index in split(seq_len(nrow(paired)), paired$cell)) {
      candidate <- paired[[paste0(metric, "_expert")]][index]
      baseline <- paired[[paste0(metric, "_reference")]][index]
      ratio <- mean(candidate) / mean(baseline)
      paired$cell_ratio[index] <- ratio
      paired$influence[index] <-
        (candidate - ratio * baseline) / mean(baseline)
    }
    family_groups <- split(
      seq_len(nrow(paired)), interaction(paired[family_names], drop = TRUE)
    )
    lapply(family_groups, function(index) {
      first_cell_rows <- index[!duplicated(paired$cell[index])]
      replication_groups <- split(index, paired$replication[index])
      replication_influence <- vapply(
        replication_groups,
        function(replication_index) mean(paired$influence[replication_index]),
        numeric(1)
      )
      data.frame(
        comparison = comparison,
        metric = metric,
        paired[index[[1L]], family_names, drop = FALSE],
        ratio = mean(paired$cell_ratio[first_cell_rows]),
        mcse = sd(replication_influence) / sqrt(length(replication_influence)),
        cells = length(first_cell_rows)
      )
    })
  })
  do.call(rbind, unlist(rows, recursive = FALSE))
}

family_ratios <- rbind(
  family_paired_ratios(pre_pr, "pre-pr"),
  family_paired_ratios(bulk, "bulk")
)
support_names <- c(
  "comparison", "metric", "support", "degree", "sample_size"
)
support_groups <- split(
  seq_len(nrow(family_ratios)),
  interaction(family_ratios[support_names], drop = TRUE)
)
support_ratios <- do.call(rbind, lapply(support_groups, function(index) {
  data.frame(
    family_ratios[index[[1L]], support_names, drop = FALSE],
    ratio = mean(family_ratios$ratio[index]),
    mcse = sqrt(sum(family_ratios$mcse[index]^2)) / length(index),
    cells = sum(family_ratios$cells[index])
  )
}))

classification <- expert[
  expert$degree == 0 &
    ((expert$support == "one-sided" & expert$scale == 1 &
        expert$direction == "lower") | expert$support == "two-sided"),
]
activation_groups <- split(
  seq_len(nrow(classification)),
  interaction(classification[c("support", "scenario", "sample_size")],
              drop = TRUE)
)
activation <- do.call(rbind, lapply(activation_groups, function(index) {
  data.frame(
    classification[index[[1L]], c("support", "scenario", "sample_size")],
    lower_rate = mean(classification$repair_lower[index]),
    upper_rate = mean(classification$repair_upper[index]),
    any_rate = mean(
      classification$repair_lower[index] | classification$repair_upper[index]
    ),
    both_rate = mean(
      classification$repair_lower[index] & classification$repair_upper[index]
    )
  )
}))

read_summaries <- function(path) {
  summary <- read.csv(path, stringsAsFactors = FALSE)
  summary$direction[is.na(summary$direction)] <- "none"
  summary
}
summary_path <- function(path) sub("\\.csv$", "-summary.csv", path)
pre_pr_summary <- read_summaries(summary_path(arguments[[1L]]))
bulk_summary <- read_summaries(summary_path(arguments[[2L]]))
expert_summary <- read_summaries(summary_path(arguments[[3L]]))
summary_keys <- keys[keys != "replication"]

component_ratios <- function(reference, comparison) {
  paired <- merge(
    expert_summary,
    reference,
    by = summary_keys,
    suffixes = c("_expert", "_reference")
  )
  metrics <- c(
    "global_squared_bias", "global_variance",
    "boundary_squared_bias", "boundary_variance"
  )
  do.call(rbind, lapply(metrics, function(metric) {
    data.frame(
      comparison = comparison,
      metric = metric,
      paired[summary_keys],
      ratio = paired[[paste0(metric, "_expert")]] /
        paired[[paste0(metric, "_reference")]]
    )
  }))
}

components <- rbind(
  component_ratios(pre_pr_summary, "pre-pr"),
  component_ratios(bulk_summary, "bulk")
)

example_targets <- data.frame(
  support = c(rep("one-sided", 4), rep("two-sided", 4)),
  scenario = c(
    "exponential", "half_normal", "gamma_075", "gamma_2",
    "uniform", "beta_1_2", "beta_075_2", "beta_2_2"
  ),
  sample_size = c(100L, 25L, 1000L, 1000L, 100L, 1000L, 1000L, 1000L),
  repair_lower = c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE),
  repair_upper = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE)
)
example_rows <- lapply(seq_len(nrow(example_targets)), function(i) {
  target <- example_targets[i, ]
  candidates <- expert[
    expert$support == target$support &
      expert$scenario == target$scenario &
      expert$sample_size == target$sample_size &
      expert$degree == 2 &
      expert$scale == 1 &
      expert$repair_lower == target$repair_lower &
      expert$repair_upper == target$repair_upper &
      expert$direction ==
        if (target$support == "one-sided") "lower" else "none",
  ]
  median_ise <- median(candidates$global_ise)
  selected <- candidates[
    which.min(abs(candidates$global_ise - median_ise)),
  ]
  data.frame(target, replication = selected$replication)
})
examples <- do.call(rbind, example_rows)

dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)
write.csv(cell_ratios, paste0(output_prefix, "-cells.csv"), row.names = FALSE)
write.csv(family_ratios, paste0(output_prefix, "-families.csv"), row.names = FALSE)
write.csv(support_ratios, paste0(output_prefix, "-support.csv"), row.names = FALSE)
write.csv(activation, paste0(output_prefix, "-activation.csv"), row.names = FALSE)
write.csv(components, paste0(output_prefix, "-components.csv"), row.names = FALSE)
write.csv(examples, paste0(output_prefix, "-examples.csv"), row.names = FALSE)
