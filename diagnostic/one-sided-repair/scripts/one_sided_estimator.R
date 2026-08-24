# Reference R implementation of the one-sided finite-boundary candidate.

one_sided_trapezoid_integral <- function(evaluation_points, values) {
  sum(diff(evaluation_points) *
    (values[-length(values)] + values[-1L]) / 2)
}

normalize_one_sided_density <- function(evaluation_points, density) {
  density <- pmax(0, density)
  density_mass <- one_sided_trapezoid_integral(evaluation_points, density)
  if (!is.finite(density_mass) || density_mass <= 0) {
    stop("density cannot be normalized")
  }
  density / density_mass
}

one_sided_smoothstep <- function(normalized_distance) {
  3 * normalized_distance^2 - 2 * normalized_distance^3
}

make_one_sided_expert_weights <- function(probabilities, sample_size,
                                          concentration = 1,
                                          shrinkage_exponent = 0.5) {
  boundary_fraction <- min(
    0.25,
    concentration / sample_size^shrinkage_exponent
  )
  boundary_weights <- 1 - one_sided_smoothstep(pmin(
    1,
    probabilities / boundary_fraction
  ))
  cbind(boundary = boundary_weights, bulk = 1 - boundary_weights)
}

one_sided_weighted_quantile <- function(values, weights, probability) {
  value_order <- order(values)
  ordered_weights <- weights[value_order]
  values[value_order][which(
    cumsum(ordered_weights) / sum(ordered_weights) >= probability
  )[[1L]]]
}

select_one_sided_bandwidth_floor <- function(
    boundary_distances,
    bandwidth_weights,
    degree,
    bandwidth_floor = c("local", "shrinking", "global", "none", "robust"),
    local_floor_fraction = 0.75) {
  bandwidth_floor <- match.arg(bandwidth_floor)
  if (bandwidth_floor == "none") {
    return(0)
  }
  if (bandwidth_floor == "global") {
    return(kde1d(boundary_distances, deg = degree)$bw)
  }
  if (bandwidth_floor == "local" || bandwidth_floor == "shrinking") {
    if (!is.finite(local_floor_fraction) || local_floor_fraction <= 0 ||
        local_floor_fraction > 1) {
      stop("local_floor_fraction must lie in (0, 1]")
    }
    local_count <- if (bandwidth_floor == "shrinking") {
      max(4L, ceiling(sqrt(length(boundary_distances))))
    } else {
      max(4L, ceiling(local_floor_fraction * length(boundary_distances)))
    }
    return(kde1d(sort(boundary_distances)[seq_len(local_count)], deg = degree)$bw)
  }
  weighted_median <- one_sided_weighted_quantile(
    boundary_distances,
    bandwidth_weights,
    0.5
  )
  weighted_mad <- one_sided_weighted_quantile(
    abs(boundary_distances - weighted_median),
    bandwidth_weights,
    0.5
  )
  if (!(weighted_mad > 0)) {
    weighted_mad <- one_sided_weighted_quantile(
      boundary_distances,
      bandwidth_weights,
      0.75
    ) - weighted_median
  }
  if (!(weighted_mad > 0)) {
    weighted_mad <- max(boundary_distances) * .Machine$double.eps
  }
  winsorized_distances <- pmin(
    boundary_distances,
    weighted_median + 6 * weighted_mad
  )
  kde1d(winsorized_distances, deg = degree)$bw
}

fit_one_sided_bulk <- function(
    observations,
    evaluation_points,
    lower_bound,
    upper_bound) {
  lower_bounded <- !is.na(lower_bound)
  bulk_model <- if (lower_bounded) {
    kde1d(observations, xmin = lower_bound)
  } else {
    kde1d(observations, xmax = upper_bound)
  }
  list(
    density = normalize_one_sided_density(
      evaluation_points,
      dkde1d(evaluation_points, bulk_model)
    ),
    boundary_probabilities = if (lower_bounded) {
      pkde1d(evaluation_points, bulk_model)
    } else {
      1 - pkde1d(evaluation_points, bulk_model)
    }
  )
}

fit_one_sided_local_polynomial <- function(
    boundary_distances,
    evaluation_distances,
    bandwidth_weights,
    minimum_bandwidth,
    degree) {
  if (!degree %in% c(1, 2)) {
    stop("local boundary degree must be 1 or 2")
  }
  selected_bandwidth <- kde1d(
    boundary_distances,
    weights = bandwidth_weights,
    deg = degree
  )$bw
  bandwidth <- max(selected_bandwidth, minimum_bandwidth)
  standardized_boundary <- evaluation_distances / bandwidth
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
  density <- numeric(length(evaluation_distances))
  for (chunk_start in seq(1L, length(evaluation_distances), by = 128L)) {
    chunk_indices <- chunk_start:min(
      chunk_start + 127L,
      length(evaluation_distances)
    )
    kernel_arguments <- outer(
      evaluation_distances[chunk_indices],
      boundary_distances,
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
    density = normalize_one_sided_density(evaluation_distances, density),
    bandwidth = bandwidth,
    selected_bandwidth = selected_bandwidth
  )
}

select_one_sided_endpoint_expert <- function(
    boundary_distances,
    tail_multiplier = 2,
    evidence_threshold = qnorm(0.95),
    minimum_finite_index = 0.9) {
  tail_count <- min(
    length(boundary_distances) - 1L,
    ceiling(tail_multiplier * sqrt(length(boundary_distances)))
  )
  ordered_distances <- sort(boundary_distances)
  tail_index <- tail_count / sum(log(
    ordered_distances[[tail_count + 1L]] /
      pmax(.Machine$double.eps, ordered_distances[seq_len(tail_count)])
  ))
  list(
    expert_choice = if (
      tail_index >= minimum_finite_index &&
        tail_index * (1 - evidence_threshold / sqrt(tail_count)) <= 1
    ) {
      "finite"
    } else {
      "bulk"
    },
    tail_index = tail_index,
    tail_count = tail_count
  )
}

fit_one_sided_candidate <- function(
    observations,
    evaluation_points,
    lower_bound = 0,
    upper_bound = NA_real_,
    concentration = 1,
    shrinkage_exponent = 0.5,
    bandwidth_floor = c("local", "shrinking", "global", "none", "robust"),
    local_floor_fraction = 0.75) {
  bandwidth_floor <- match.arg(bandwidth_floor)
  if (is.na(lower_bound) == is.na(upper_bound)) {
    stop("exactly one finite boundary must be supplied")
  }
  if (length(observations) < 4L || any(!is.finite(observations)) ||
      any(!is.finite(evaluation_points)) ||
      is.unsorted(evaluation_points, strictly = TRUE)) {
    stop("observations and evaluation_points must be finite and ordered")
  }
  lower_bounded <- !is.na(lower_bound)
  boundary_distances <- if (lower_bounded) {
    observations - lower_bound
  } else {
    upper_bound - observations
  }
  if (any(boundary_distances < 0)) {
    stop("observations must lie inside the supplied support")
  }
  evaluation_distances <- if (lower_bounded) {
    evaluation_points - lower_bound
  } else {
    rev(upper_bound - evaluation_points)
  }
  bulk_fit <- fit_one_sided_bulk(
    observations,
    evaluation_points,
    lower_bound,
    upper_bound
  )
  selection <- select_one_sided_endpoint_expert(boundary_distances)
  if (selection$expert_choice == "bulk") {
    return(list(
      density = bulk_fit$density,
      expert_choice = "bulk",
      tail_index = selection$tail_index,
      tail_count = selection$tail_count,
      boundary_bandwidths = c(linear = NA_real_, quadratic = NA_real_),
      selected_boundary_bandwidths = c(
        linear = NA_real_, quadratic = NA_real_
      ),
      bandwidth_floor = bandwidth_floor,
      local_floor_fraction = local_floor_fraction,
      boundary_fraction = min(
        0.25,
        concentration / length(observations)^shrinkage_exponent
      )
    ))
  }
  probability_ranks <- (rank(
    boundary_distances,
    ties.method = "average"
  ) - 0.5) / length(observations)
  bandwidth_weights <- make_one_sided_expert_weights(
    probability_ranks,
    length(observations),
    concentration,
    shrinkage_exponent
  )[, "boundary"]
  local_fits <- lapply(c(1L, 2L), function(degree) {
    fit_one_sided_local_polynomial(
      boundary_distances,
      evaluation_distances,
      bandwidth_weights,
      select_one_sided_bandwidth_floor(
        boundary_distances,
        bandwidth_weights,
        degree,
        bandwidth_floor,
        local_floor_fraction
      ),
      degree
    )
  })
  boundary_density <- rowMeans(vapply(
    local_fits,
    function(local_fit) local_fit$density,
    numeric(length(evaluation_distances))
  ))
  if (!lower_bounded) {
    boundary_density <- rev(boundary_density)
  }
  evaluation_weights <- make_one_sided_expert_weights(
    bulk_fit$boundary_probabilities,
    length(observations),
    concentration,
    shrinkage_exponent
  )
  list(
    density = normalize_one_sided_density(
      evaluation_points,
      evaluation_weights[, "boundary"] * boundary_density +
        evaluation_weights[, "bulk"] * bulk_fit$density
    ),
    expert_choice = "finite",
    tail_index = selection$tail_index,
    tail_count = selection$tail_count,
    boundary_bandwidths = setNames(
      vapply(local_fits, function(local_fit) local_fit$bandwidth, numeric(1L)),
      c("linear", "quadratic")
    ),
    selected_boundary_bandwidths = setNames(
      vapply(
        local_fits,
        function(local_fit) local_fit$selected_bandwidth,
        numeric(1L)
      ),
      c("linear", "quadratic")
    ),
    bandwidth_floor = bandwidth_floor,
    local_floor_fraction = local_floor_fraction,
    boundary_fraction = min(
      0.25,
      concentration / length(observations)^shrinkage_exponent
    )
  )
}
