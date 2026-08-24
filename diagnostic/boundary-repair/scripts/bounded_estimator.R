# Reference R implementation of the selected two-sided boundary estimator.

trapezoid_integral <- function(evaluation_points, values) {
  sum(diff(evaluation_points) *
    (values[-length(values)] + values[-1L]) / 2)
}

normalize_density <- function(evaluation_points, density) {
  density <- pmax(0, density)
  density_mass <- trapezoid_integral(evaluation_points, density)
  if (!is.finite(density_mass) || density_mass <= 0) {
    stop("density cannot be normalized")
  }
  density / density_mass
}

smoothstep <- function(normalized_distance) {
  3 * normalized_distance^2 - 2 * normalized_distance^3
}

make_bounded_expert_weights <- function(probabilities, sample_size,
                                        concentration = 1,
                                        shrinkage_exponent = 0.5) {
  boundary_fraction <- min(
    0.25,
    concentration / sample_size^shrinkage_exponent
  )
  lower_weights <- 1 - smoothstep(pmin(1, probabilities / boundary_fraction))
  upper_weights <- 1 - smoothstep(pmin(
    1,
    (1 - probabilities) / boundary_fraction
  ))
  cbind(
    lower = lower_weights,
    bulk = pmax(0, 1 - lower_weights - upper_weights),
    upper = upper_weights
  )
}

fit_local_polynomial_boundary <- function(
    observations,
    evaluation_points,
    bandwidth_weights,
    minimum_bandwidth,
    degree) {
  if (!degree %in% c(1, 2)) {
    stop("local boundary degree must be 1 or 2")
  }
  selected_bandwidth <- kde1d(
    observations,
    weights = bandwidth_weights,
    deg = degree
  )$bw
  bandwidth <- max(selected_bandwidth, minimum_bandwidth)
  standardized_boundary <- evaluation_points / bandwidth
  truncated_moment_0 <- pnorm(standardized_boundary)
  truncated_moment_1 <- -dnorm(standardized_boundary)
  truncated_moment_2 <- truncated_moment_0 -
    standardized_boundary * dnorm(standardized_boundary)
  if (degree == 1) {
    moment_determinant <- truncated_moment_0 * truncated_moment_2 -
      truncated_moment_1^2
    equivalent_kernel_coefficients <- cbind(
      truncated_moment_2 / moment_determinant,
      -truncated_moment_1 / moment_determinant
    )
  } else {
    truncated_moment_3 <- -(
      standardized_boundary^2 + 2
    ) * dnorm(standardized_boundary)
    truncated_moment_4 <- 3 * truncated_moment_0 - (
      standardized_boundary^3 + 3 * standardized_boundary
    ) * dnorm(standardized_boundary)
    moment_determinant <-
      truncated_moment_0 * (
        truncated_moment_2 * truncated_moment_4 - truncated_moment_3^2
      ) -
      truncated_moment_1 * (
        truncated_moment_1 * truncated_moment_4 -
          truncated_moment_2 * truncated_moment_3
      ) +
      truncated_moment_2 * (
        truncated_moment_1 * truncated_moment_3 - truncated_moment_2^2
      )
    equivalent_kernel_coefficients <- cbind(
      (
        truncated_moment_2 * truncated_moment_4 - truncated_moment_3^2
      ) / moment_determinant,
      (
        truncated_moment_2 * truncated_moment_3 -
          truncated_moment_1 * truncated_moment_4
      ) / moment_determinant,
      (
        truncated_moment_1 * truncated_moment_3 - truncated_moment_2^2
      ) / moment_determinant
    )
  }
  density <- numeric(length(evaluation_points))
  for (chunk_start in seq(1L, length(evaluation_points), by = 128L)) {
    chunk_indices <- chunk_start:min(
      chunk_start + 127L,
      length(evaluation_points)
    )
    kernel_arguments <- outer(
      evaluation_points[chunk_indices],
      observations,
      "-"
    ) / bandwidth
    base_kernels <- dnorm(kernel_arguments)
    corrected_kernels <- base_kernels *
      equivalent_kernel_coefficients[chunk_indices, 1L]
    for (polynomial_order in seq_len(degree)) {
      corrected_kernels <- corrected_kernels +
        base_kernels *
          equivalent_kernel_coefficients[
            chunk_indices,
            polynomial_order + 1L
          ] * kernel_arguments^polynomial_order
    }
    density[chunk_indices] <- rowMeans(corrected_kernels) / bandwidth
  }
  list(
    density = normalize_density(evaluation_points, density),
    bandwidth = bandwidth,
    selected_bandwidth = selected_bandwidth
  )
}

fit_bounded_finite_expert <- function(observations, evaluation_points,
                                      concentration = 1,
                                      shrinkage_exponent = 0.5) {
  probability_ranks <- (rank(observations, ties.method = "average") - 0.5) /
    length(observations)
  bandwidth_weights <- make_bounded_expert_weights(
    probability_ranks,
    length(observations),
    concentration,
    shrinkage_exponent
  )
  bulk_model <- kde1d(observations, xmin = 0, xmax = 1)
  bulk_density <- normalize_density(
    evaluation_points,
    dkde1d(evaluation_points, bulk_model)
  )
  endpoint_densities <- lapply(c(1L, 2L), function(degree) {
    minimum_bandwidth <- kde1d(observations, deg = degree)$bw
    lower_density <- fit_local_polynomial_boundary(
      observations,
      evaluation_points,
      bandwidth_weights[, "lower"],
      minimum_bandwidth,
      degree
    )$density
    upper_density <- rev(fit_local_polynomial_boundary(
      1 - observations,
      rev(1 - evaluation_points),
      bandwidth_weights[, "upper"],
      minimum_bandwidth,
      degree
    )$density)
    cbind(lower = lower_density, bulk = bulk_density, upper = upper_density)
  })
  component_densities <- (
    endpoint_densities[[1L]] + endpoint_densities[[2L]]
  ) / 2
  evaluation_weights <- make_bounded_expert_weights(
    pkde1d(evaluation_points, bulk_model),
    length(observations),
    concentration,
    shrinkage_exponent
  )
  list(
    density = normalize_density(
      evaluation_points,
      rowSums(evaluation_weights * component_densities)
    ),
    component_densities = component_densities,
    evaluation_weights = evaluation_weights,
    boundary_fraction = min(
      0.25,
      concentration / length(observations)^shrinkage_exponent
    )
  )
}

estimate_bounded_tail_indices <- function(observations, tail_counts) {
  if (any(tail_counts < 1L) || any(tail_counts >= length(observations))) {
    stop("tail_counts must be between one and sample size minus one")
  }
  ordered_observations <- sort(observations)
  upper_distances <- sort(1 - observations)
  tail_indices <- vapply(tail_counts, function(tail_count) {
    c(
      lower = tail_count / sum(log(
        ordered_observations[[tail_count + 1L]] /
          pmax(
            .Machine$double.eps,
            ordered_observations[seq_len(tail_count)]
          )
      )),
      upper = tail_count / sum(log(
        upper_distances[[tail_count + 1L]] /
          pmax(.Machine$double.eps, upper_distances[seq_len(tail_count)])
      ))
    )
  }, numeric(2L))
  if (length(tail_counts) == 1L) {
    tail_indices <- matrix(
      tail_indices,
      nrow = 2L,
      dimnames = list(c("lower", "upper"), NULL)
    )
  }
  rownames(tail_indices) <- c("lower", "upper")
  colnames(tail_indices) <- as.character(tail_counts)
  tail_indices
}

select_bounded_endpoint_experts_single_scale <- function(
    observations,
    tail_multiplier = 1,
    evidence_threshold = qnorm(0.95),
    minimum_finite_index = 0.9) {
  if (!is.finite(tail_multiplier) || tail_multiplier <= 0) {
    stop("tail_multiplier must be finite and positive")
  }
  tail_count <- min(
    length(observations) - 1L,
    ceiling(tail_multiplier * sqrt(length(observations)))
  )
  tail_indices <- estimate_bounded_tail_indices(
    observations,
    tail_count
  )[, 1L]
  expert_choices <- rep("bulk", length(tail_indices))
  names(expert_choices) <- names(tail_indices)
  expert_choices[
    tail_indices >= minimum_finite_index &
      tail_indices * (1 - evidence_threshold / sqrt(tail_count)) <= 1
  ] <- "finite"
  list(
    expert_choices = expert_choices,
    tail_indices = tail_indices,
    tail_count = tail_count,
    evidence_threshold = evidence_threshold,
    minimum_finite_index = minimum_finite_index
  )
}

select_bounded_endpoint_experts_multiscale <- function(
    observations,
    evidence_threshold = qnorm(0.90),
    tail_multipliers = c(1, 1.5, 2),
    maximum_tail_fraction = 0.15,
    rule = c("vote", "median")) {
  rule <- match.arg(rule)
  tail_counts <- pmin(
    floor(maximum_tail_fraction * length(observations)),
    ceiling(tail_multipliers * sqrt(length(observations)))
  )
  tail_counts <- pmax(1L, pmin(length(observations) - 1L, tail_counts))
  tail_indices <- estimate_bounded_tail_indices(observations, tail_counts)
  standardized_indices <- sweep(
    log(tail_indices),
    2L,
    sqrt(tail_counts),
    "*"
  )
  expert_choices <- apply(standardized_indices, 1L, function(scores) {
    if (rule == "median") {
      median_score <- median(scores)
      if (abs(median_score) <= evidence_threshold) {
        "finite"
      } else {
        "bulk"
      }
    } else {
      required_votes <- floor(length(scores) / 2) + 1L
      if (sum(abs(scores) <= evidence_threshold) >= required_votes) {
        "finite"
      } else {
        "bulk"
      }
    }
  })
  list(
    expert_choices = expert_choices,
    tail_indices = tail_indices,
    standardized_indices = standardized_indices,
    tail_counts = tail_counts,
    evidence_threshold = evidence_threshold,
    tail_multipliers = tail_multipliers,
    maximum_tail_fraction = maximum_tail_fraction,
    rule = rule
  )
}

select_bounded_endpoint_experts <- function(
    observations,
    evidence_threshold = qnorm(0.95),
    minimum_finite_index = 0.9) {
  select_bounded_endpoint_experts_single_scale(
    observations,
    tail_multiplier = 2,
    evidence_threshold,
    minimum_finite_index
  )
}

fit_bounded_estimator_components <- function(
    observations,
    evaluation_points,
    concentration = 1,
    shrinkage_exponent = 0.5) {
  if (length(observations) < 4L || any(!is.finite(observations)) ||
      any(observations < 0 | observations > 1)) {
    stop("observations must contain at least four finite values in [0, 1]")
  }
  if (any(!is.finite(evaluation_points)) ||
      is.unsorted(evaluation_points, strictly = TRUE) ||
      any(evaluation_points < 0 | evaluation_points > 1)) {
    stop("evaluation_points must be strictly increasing in [0, 1]")
  }
  finite_fit <- fit_bounded_finite_expert(
    observations,
    evaluation_points,
    concentration,
    shrinkage_exponent
  )
  selection <- select_bounded_endpoint_experts(observations)
  list(
    finite_fit = finite_fit,
    selection = selection
  )
}

fuse_bounded_estimator_components <- function(
    evaluation_points,
    estimator_components) {
  finite_fit <- estimator_components$finite_fit
  selection <- estimator_components$selection
  component_densities <- finite_fit$component_densities
  for (boundary in c("lower", "upper")) {
    if (selection$expert_choices[[boundary]] != "finite") {
      component_densities[, boundary] <- component_densities[, "bulk"]
    }
  }
  list(
    density = normalize_density(
      evaluation_points,
      rowSums(finite_fit$evaluation_weights * component_densities)
    ),
    component_densities = component_densities,
    evaluation_weights = finite_fit$evaluation_weights,
    expert_choices = selection$expert_choices,
    tail_indices = selection$tail_indices,
    tail_count = selection$tail_count,
    boundary_fraction = finite_fit$boundary_fraction
  )
}

fit_bounded_selected_estimator <- function(
    observations,
    evaluation_points,
    concentration = 1,
    shrinkage_exponent = 0.5) {
  fuse_bounded_estimator_components(
    evaluation_points,
    fit_bounded_estimator_components(
      observations,
      evaluation_points,
      concentration,
      shrinkage_exponent
    )
  )
}
