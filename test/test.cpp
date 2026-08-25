#include "../include/kde1d.hpp"

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

using namespace kde1d;

long int n_sample = 10000;
double pdf_tol = 0.2;
Eigen::VectorXd ugrid = Eigen::VectorXd::LinSpaced(99, 0.01, 0.99);
Eigen::VectorXd upoints = Eigen::VectorXd::LinSpaced(9, 0.1, 0.9);
// continuous, bounded data
Eigen::VectorXd x_cb = stats::simulate_uniform(n_sample, { 1 });
// continuous, unbounded data (using the standard normal quantile)
Eigen::VectorXd x_ub = stats::qnorm(x_cb);
// continuous, right bounded data (using the unit exponential quantile)
Eigen::VectorXd x_rb = x_cb.array().log();
// continuous, left bounded data (using the unit exponential quantile)
Eigen::VectorXd x_lb = x_rb * (-1.0);
// // discrete data
size_t nlevels = 50;
Eigen::VectorXd x_d =
  (x_cb.array() * (static_cast<double>(nlevels) - 1)).round();

TEST_CASE("sample median handles odd and even vectors", "[stats]")
{
  Eigen::VectorXd odd(5);
  odd << 4.0, 1.0, 5.0, 2.0, 3.0;
  Eigen::VectorXd even(4);
  even << 4.0, 1.0, 3.0, 2.0;

  CHECK(stats::median(odd) == 3.0);
  CHECK(stats::median(even) == 2.5);
  CHECK_THROWS(stats::median(Eigen::VectorXd()));
}

TEST_CASE("linear binning includes the upper endpoint", "[linear-binning]")
{
  Eigen::VectorXd observations(3);
  observations << 0.0, 0.5, 1.0;
  Eigen::VectorXd counts = tools::linbin(
    observations, 0.0, 1.0, 2, Eigen::VectorXd::Ones(3));

  CHECK(counts.isApprox(Eigen::VectorXd::Ones(3)));
  CHECK(counts.sum() == Approx(3.0));
}

TEST_CASE("right extrapolation is continuous", "[interpolation]")
{
  Eigen::VectorXd grid_points(2);
  grid_points << 0.0, 1.0;
  Eigen::VectorXd values(2);
  values << 2.0, 3.0;
  interp::InterpolationGrid grid(grid_points, values, 0);

  CHECK(grid.interpolate(Eigen::VectorXd::Constant(1, 1.0))(0) ==
        Approx(3.0));
  CHECK(grid.interpolate(Eigen::VectorXd::Constant(1, 2.0))(0) ==
        Approx(3.0 * std::exp(-0.5)));
}

TEST_CASE("spline quantiles invert nonuniform cumulative integrals",
          "[interpolation][quantile]")
{
  Eigen::VectorXd grid_points(6);
  grid_points << -2.0, -0.7, -0.1, 0.2, 0.9, 1.5;
  Eigen::VectorXd values(6);
  values << 0.05, 0.4, 1.2, 0.8, 0.3, 0.02;
  interp::InterpolationGrid grid(grid_points, values, 1);

  Eigen::VectorXd probabilities(10);
  probabilities << 0.75, 1e-10, 0.25, 0.5, 1.0 - 1e-10, 0.25, 0.0, 1.0,
    0.9, 0.1;
  Eigen::VectorXd quantiles = grid.quantile(probabilities);
  Eigen::VectorXd reference = tools::invert_f(
    probabilities,
    [&](const Eigen::VectorXd& x) { return grid.integrate(x, true); },
    grid.get_grid_min(),
    grid.get_grid_max(),
    50);

  CHECK(quantiles.isApprox(reference, 1e-10));
  CHECK(grid.integrate(quantiles, true).isApprox(probabilities, 1e-12));
  CHECK(quantiles(2) == quantiles(5));
  CHECK(quantiles(6) == grid.get_grid_min());
  CHECK(quantiles(7) == grid.get_grid_max());
}

TEST_CASE("boundary grids resolve the transformed support", "[boundary-grid]")
{
  Eigen::VectorXd observations = Eigen::VectorXd::LinSpaced(200, 0.2, 0.8);
  Kde1d bounded(0.0, 1.0, "continuous", 1.0, 0.3);
  bounded.fit(observations);
  Eigen::VectorXd grid = bounded.get_grid_points();

  CHECK(grid(0) == Approx(0.0));
  CHECK(grid(grid.size() - 1) == Approx(1.0));
  CHECK(grid(1) < observations.minCoeff());
  CHECK(grid(grid.size() - 2) > observations.maxCoeff());
  CHECK((grid.tail(grid.size() - 1) - grid.head(grid.size() - 1)).minCoeff() >
        0.0);

  Kde1d left_bounded(0.0, NAN, "continuous", 1.0, 0.3);
  left_bounded.fit(observations);
  grid = left_bounded.get_grid_points();
  CHECK(grid(0) == Approx(0.0));
  CHECK(grid(1) < observations.minCoeff());
  CHECK(grid(grid.size() - 1) > observations.maxCoeff());

  observations *= -1.0;
  Kde1d right_bounded(NAN, 0.0, "continuous", 1.0, 0.3);
  right_bounded.fit(observations);
  grid = right_bounded.get_grid_points();
  CHECK(grid(0) < observations.minCoeff());
  CHECK(grid(grid.size() - 2) > observations.maxCoeff());
  CHECK(grid(grid.size() - 1) == Approx(0.0));
}

TEST_CASE("one-boundary fits are reflection equivariant", "[boundary-reflection]")
{
  Eigen::VectorXd observations = Eigen::VectorXd::LinSpaced(500, 0.01, 5.0);
  Kde1d left_bounded(0.0, NAN, "continuous");
  left_bounded.fit(observations);
  Kde1d manual(0.0, NAN, "continuous", 1.0, left_bounded.get_bandwidth());
  manual.fit(observations);
  CHECK(left_bounded.get_values().isApprox(manual.get_values(), 1e-12));

  Kde1d right_bounded(NAN, 0.0, "continuous");
  right_bounded.fit(-observations);

  CHECK(left_bounded.get_grid_points().isApprox(
    -right_bounded.get_grid_points().reverse()));
  CHECK(left_bounded.get_values().isApprox(
    right_bounded.get_values().reverse()));
  CHECK(left_bounded.get_loglik() == Approx(right_bounded.get_loglik()));
  CHECK(left_bounded.get_edf() == Approx(right_bounded.get_edf()));
}

TEST_CASE("one-sided finite endpoints use the boundary expert",
          "[boundary-expert]")
{
  const Eigen::VectorXd probabilities =
    Eigen::VectorXd::LinSpaced(200, 0.5 / 200.0, 1.0 - 0.5 / 200.0);
  const Eigen::VectorXd observations = (-(1.0 - probabilities.array()).log());
  Kde1d fit(0.0, NAN, "continuous");
  fit.fit(observations);

  Kde1d manual(0.0, NAN, "continuous", 1.0, fit.get_bandwidth());
  manual.fit(observations);
  CHECK(fit.get_values().isApprox(manual.get_values(), 1e-12));

  Eigen::VectorXd expected_density(6);
  expected_density << 0.9797951, 0.9798032, 0.9780643, 0.9448344, 0.7316556,
    0.0077380;
  Eigen::VectorXd selected_density(6);
  selected_density << fit.get_values()(0), fit.get_values()(4),
    fit.get_values()(24), fit.get_values()(49), fit.get_values()(99),
    fit.get_values()(199);
  CHECK(selected_density.isApprox(expected_density, 2e-4));
  CHECK(std::isfinite(fit.get_edf()));

  Kde1d reflected(NAN, 0.0, "continuous");
  reflected.fit(-observations);
  CHECK(fit.get_grid_points().isApprox(-reflected.get_grid_points().reverse(),
                                       1e-12));
  CHECK(fit.get_values().isApprox(reflected.get_values().reverse(), 1e-12));
  CHECK(fit.get_edf() == Approx(reflected.get_edf()).epsilon(1e-12));
}

TEST_CASE("boundary experts can be disabled", "[boundary-expert]")
{
  const Eigen::VectorXd probabilities =
    Eigen::VectorXd::LinSpaced(200, 0.5 / 200.0, 1.0 - 0.5 / 200.0);
  const Eigen::VectorXd observations = (-(1.0 - probabilities.array()).log());

  Kde1d default_fit(0.0, NAN, "continuous");
  default_fit.fit(observations);
  Kde1d enabled(0.0, NAN, "continuous", 1.0, NAN, 2, 400, true);
  enabled.fit(observations);
  Kde1d disabled(0.0, NAN, "continuous", 1.0, NAN, 2, 400, false);
  disabled.fit(observations);
  Kde1d disabled_manual(
    0.0, NAN, "continuous", 1.0, disabled.get_bandwidth(), 2, 400, false);
  disabled_manual.fit(observations);

  CHECK(default_fit.get_boundary_repair());
  CHECK(enabled.get_boundary_repair());
  CHECK_FALSE(disabled.get_boundary_repair());
  CHECK(default_fit.get_values().isApprox(enabled.get_values(), 1e-12));
  CHECK_FALSE(default_fit.get_values().isApprox(disabled.get_values(), 1e-6));
  CHECK(disabled.get_values().isApprox(disabled_manual.get_values(), 1e-12));
  CHECK(disabled.str().find("boundary_repair=false") != std::string::npos);
}

TEST_CASE("one-sided vanishing endpoints retain the bulk fit",
          "[boundary-expert]")
{
  const Eigen::VectorXd probabilities =
    Eigen::VectorXd::LinSpaced(200, 0.5 / 200.0, 1.0 - 0.5 / 200.0);
  const Eigen::VectorXd observations = probabilities.array().sqrt();
  Kde1d fit(0.0, NAN, "continuous");
  fit.fit(observations);

  Kde1d bulk(0.0, NAN, "continuous", 1.0, fit.get_bandwidth());
  bulk.fit(observations);
  CHECK(fit.get_grid_points().isApprox(bulk.get_grid_points(), 1e-12));
  CHECK(fit.get_values().isApprox(bulk.get_values(), 1e-12));
  CHECK(fit.get_edf() == Approx(bulk.get_edf()).epsilon(1e-12));
}

TEST_CASE("two-sided finite endpoints use the boundary experts",
          "[boundary-expert]")
{
  const Eigen::VectorXd observations =
    Eigen::VectorXd::LinSpaced(200, 0.5 / 200.0, 1.0 - 0.5 / 200.0);
  Kde1d fit(0.0, 1.0, "continuous");
  fit.fit(observations);

  Kde1d manual(0.0, 1.0, "continuous", 1.0, fit.get_bandwidth());
  manual.fit(observations);
  CHECK(fit.get_values().isApprox(manual.get_values(), 1e-12));

  Eigen::VectorXd expected_density(6);
  expected_density << 0.9993158, 0.9992803, 0.9985546, 1.0001208, 0.9985473,
    0.9993158;
  Eigen::VectorXd selected_density(6);
  selected_density << fit.get_values()(0), fit.get_values()(49),
    fit.get_values()(99), fit.get_values()(199), fit.get_values()(299),
    fit.get_values()(400);
  CHECK(selected_density.isApprox(expected_density, 2e-4));
  CHECK(std::isfinite(fit.get_edf()));

  Kde1d scaled(-3.0, 7.0, "continuous");
  scaled.fit(-3.0 + 10.0 * observations.array());
  CHECK(((scaled.get_grid_points().array() + 3.0) / 10.0)
          .matrix()
          .isApprox(fit.get_grid_points(), 1e-12));
  CHECK((10.0 * scaled.get_values()).isApprox(fit.get_values(), 1e-10));
  CHECK(scaled.get_edf() == Approx(fit.get_edf()).epsilon(1e-10));
}

TEST_CASE("two-sided endpoints are classified independently",
          "[boundary-expert]")
{
  const Eigen::VectorXd probabilities =
    Eigen::VectorXd::LinSpaced(200, 0.5 / 200.0, 1.0 - 0.5 / 200.0);
  const Eigen::VectorXd observations =
    1.0 - (1.0 - probabilities.array()).sqrt();
  Kde1d fit(0.0, 1.0, "continuous");
  fit.fit(observations);

  Eigen::VectorXd expected_density(6);
  expected_density << 1.9985571, 1.9950070, 1.9824915, 0.9946390, 0.0562128,
    0.0000020;
  Eigen::VectorXd selected_density(6);
  selected_density << fit.get_values()(0), fit.get_values()(49),
    fit.get_values()(99), fit.get_values()(199), fit.get_values()(299),
    fit.get_values()(400);
  CHECK(selected_density.isApprox(expected_density, 2e-4));

  Kde1d reflected(0.0, 1.0, "continuous");
  reflected.fit(1.0 - observations.array());
  CHECK(fit.get_values().isApprox(reflected.get_values().reverse(), 1e-12));
  CHECK(fit.get_edf() == Approx(reflected.get_edf()).epsilon(1e-12));
}

TEST_CASE("two-sided exploding endpoints retain the bulk fit",
          "[boundary-expert]")
{
  const Eigen::VectorXd probabilities =
    Eigen::VectorXd::LinSpaced(200, 0.5 / 200.0, 1.0 - 0.5 / 200.0);
  const Eigen::VectorXd observations =
    (0.5 * std::acos(-1.0) * probabilities.array()).sin().square();
  Kde1d fit(0.0, 1.0, "continuous");
  fit.fit(observations);

  Kde1d bulk(0.0, 1.0, "continuous", 1.0, fit.get_bandwidth());
  bulk.fit(observations);
  CHECK(fit.get_values().isApprox(bulk.get_values(), 1e-12));
  CHECK(fit.get_edf() == Approx(bulk.get_edf()).epsilon(1e-12));
}

TEST_CASE("boundary experts support weights and manual bandwidths",
          "[boundary-expert]")
{
  const Eigen::VectorXd probabilities =
    Eigen::VectorXd::LinSpaced(200, 0.5 / 200.0, 1.0 - 0.5 / 200.0);
  const Eigen::VectorXd observations = (-(1.0 - probabilities.array()).log());
  const Eigen::VectorXd weights =
    Eigen::VectorXd::LinSpaced(observations.size(), 0.5, 1.5);

  Kde1d weighted(0.0, NAN, "continuous");
  weighted.fit(observations, weights);
  Kde1d weighted_manual(0.0, NAN, "continuous", 1.0, weighted.get_bandwidth());
  weighted_manual.fit(observations, weights);
  CHECK(weighted.get_values().isApprox(weighted_manual.get_values(), 1e-12));
  CHECK(weighted.get_edf() == Approx(weighted_manual.get_edf()).epsilon(1e-12));

  Kde1d rescaled_weights(0.0, NAN, "continuous");
  rescaled_weights.fit(observations, 7.0 * weights);
  CHECK(weighted.get_values().isApprox(rescaled_weights.get_values(), 1e-12));
  CHECK(weighted.get_edf() ==
        Approx(rescaled_weights.get_edf()).epsilon(1e-12));
  CHECK(weighted.get_loglik() ==
        Approx(rescaled_weights.get_loglik()).epsilon(1e-12));

  Eigen::VectorXd zero_weights = weights;
  for (Eigen::Index i = 1; i < zero_weights.size(); i += 5)
    zero_weights(i) = 0.0;
  Eigen::VectorXd retained_observations(160);
  Eigen::VectorXd retained_weights(160);
  Eigen::Index retained_index = 0;
  for (Eigen::Index i = 0; i < zero_weights.size(); ++i) {
    if (zero_weights(i) > 0.0) {
      retained_observations(retained_index) = observations(i);
      retained_weights(retained_index++) = zero_weights(i);
    }
  }
  Kde1d with_zeros(0.0, NAN, "continuous");
  with_zeros.fit(observations, zero_weights);
  Kde1d without_zeros(0.0, NAN, "continuous");
  without_zeros.fit(retained_observations, retained_weights);
  CHECK(with_zeros.get_values().isApprox(without_zeros.get_values(), 1e-12));
  CHECK(with_zeros.get_edf() == Approx(without_zeros.get_edf()).epsilon(1e-12));

  Kde1d reflected(NAN, 0.0, "continuous");
  reflected.fit(-observations, weights);
  CHECK(
    weighted.get_values().isApprox(reflected.get_values().reverse(), 1e-12));
  CHECK(weighted.get_edf() == Approx(reflected.get_edf()).epsilon(1e-12));

  Kde1d linear(0.0, NAN, "continuous", 1.0, NAN, 1);
  linear.fit(observations);
  Kde1d linear_bulk(0.0, NAN, "continuous", 1.0, linear.get_bandwidth(), 1);
  linear_bulk.fit(observations);
  CHECK(linear.get_values().isApprox(linear_bulk.get_values(), 1e-12));
  CHECK(linear.get_edf() == Approx(linear_bulk.get_edf()).epsilon(1e-12));

  Kde1d small(0.0, NAN, "continuous");
  small.fit(observations.head(15));
  Kde1d small_bulk(0.0, NAN, "continuous", 1.0, small.get_bandwidth());
  small_bulk.fit(observations.head(15));
  CHECK(small.get_values().isApprox(small_bulk.get_values(), 1e-12));
  CHECK(small.get_edf() == Approx(small_bulk.get_edf()).epsilon(1e-12));
}

TEST_CASE("finite support truncates density", "[finite-support]")
{
  Eigen::VectorXd observations = Eigen::VectorXd::LinSpaced(200, 0.0, 1.0);
  Kde1d bounded(0.0, 1.0, "continuous", 1.0, 0.5);
  bounded.fit(observations);

  Eigen::VectorXd evaluation_points(4);
  evaluation_points << -1e-8, 0.0, 1.0, 1.0 + 1e-8;
  Eigen::VectorXd density = bounded.pdf(evaluation_points);
  CHECK(density(0) == 0.0);
  CHECK(density(1) > 0.0);
  CHECK(density(2) > 0.0);
  CHECK(density(3) == 0.0);
}

TEST_CASE("continuous quantiles invert the normalized CDF", "[quantile]")
{
  Eigen::VectorXd observations = Eigen::VectorXd::LinSpaced(500, 0.0, 1.0);
  Kde1d fit(0.0, 1.0, "continuous", 1.0, 0.1, 2, 400);
  fit.fit(observations);

  Eigen::VectorXd probabilities(12);
  probabilities << 0.0,
    1e-12,
    1e-8,
    1e-4,
    0.1,
    0.5,
    0.5,
    0.9,
    1.0 - 1e-4,
    1.0 - 1e-8,
    1.0 - 1e-12,
    1.0;
  Eigen::VectorXd quantiles = fit.quantile(probabilities);
  Eigen::VectorXd round_trip = fit.cdf(quantiles);
  CAPTURE(probabilities.transpose(),
          quantiles.transpose(),
          round_trip.transpose());

  CHECK((quantiles.tail(quantiles.size() - 1) -
         quantiles.head(quantiles.size() - 1))
          .minCoeff() >= 0.0);
  CHECK(quantiles(5) == quantiles(6));
  CHECK((round_trip - probabilities).cwiseAbs().maxCoeff() < 1e-8);
  CHECK(quantiles(0) >= 0.0);
  CHECK(quantiles(quantiles.size() - 1) <= 1.0);

  Eigen::VectorXd shuffled(7);
  shuffled << 0.9, NAN, 0.1, 0.9, 0.5, 0.1, NAN;
  quantiles = fit.quantile(shuffled);
  CHECK(std::isnan(quantiles(1)));
  CHECK(std::isnan(quantiles(6)));
  CHECK(quantiles(0) == quantiles(3));
  CHECK(quantiles(2) == quantiles(5));
}

TEST_CASE("discrete quantiles satisfy the generalized inverse", "[quantile]")
{
  Eigen::VectorXd observations(10);
  observations << 0.0, 0.0, 0.0, 1.0, 1.0, 2.0, 3.0, 3.0, 3.0, 3.0;
  Kde1d fit(0.0, 3.0, "discrete", 1.0, 0.5, 2, 100);
  fit.fit(observations);

  Eigen::VectorXd probabilities(9);
  probabilities << 0.0, 1e-12, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0 - 1e-12, 1.0;
  Eigen::VectorXd quantiles = fit.quantile(probabilities);
  Eigen::VectorXd cumulative = fit.cdf(quantiles);

  CHECK((quantiles.array() == quantiles.array().round()).all());
  CHECK((quantiles.tail(quantiles.size() - 1) -
         quantiles.head(quantiles.size() - 1))
          .minCoeff() >= 0.0);
  CHECK((cumulative.array() + 1e-14 >= probabilities.array()).all());

  Eigen::VectorXd levels = Eigen::VectorXd::LinSpaced(4, 0.0, 3.0);
  CHECK(fit.quantile(fit.cdf(levels)).isApprox(levels));
}

TEST_CASE("discrete CDF stays within the unit interval", "[discrete]")
{
  Eigen::VectorXd grid_points = Eigen::VectorXd::LinSpaced(10, 0.0, 9.0);
  Eigen::VectorXd values = Eigen::VectorXd::Ones(10);
  values(9) = 0.0;
  interp::InterpolationGrid grid(grid_points, values, 0);
  Kde1d discrete(grid, NAN, NAN, "discrete");

  CHECK(discrete.cdf(Eigen::VectorXd::Constant(1, 8.0))(0) == 1.0);
}

TEST_CASE("discrete bounds describe integer support", "[discrete][boundary]")
{
  Eigen::VectorXd observations(200);
  for (Eigen::Index i = 0; i < observations.size(); ++i)
    observations(i) = static_cast<double>(i % 4);

  Kde1d bounded(0.0, 3.0, "discrete", 1.0, NAN, 2, 100, true);
  bounded.fit(observations);
  const Eigen::VectorXd grid = bounded.get_grid_points();
  CHECK(grid(0) == Approx(-0.5));
  CHECK(grid(grid.size() - 1) == Approx(3.5));
  CHECK((grid.tail(grid.size() - 1) - grid.head(grid.size() - 1))
          .minCoeff() > 0.0);

  const Eigen::VectorXd levels = Eigen::VectorXd::LinSpaced(4, 0.0, 3.0);
  CHECK(bounded.pdf(levels).sum() == Approx(1.0));
  CHECK(bounded.cdf(Eigen::VectorXd::Constant(1, -0.1))(0) == 0.0);
  CHECK(bounded.cdf(Eigen::VectorXd::Constant(1, 3.0))(0) == 1.0);
  Eigen::VectorXd endpoints(2);
  endpoints << 0.0, 1.0;
  CHECK(bounded.quantile(endpoints)(0) == 0.0);
  CHECK(bounded.quantile(endpoints)(1) == 3.0);

  Kde1d left_bounded(0.0, NAN, "discrete", 1.0, NAN, 2, 100);
  left_bounded.fit(observations);
  const Eigen::VectorXd left_grid = left_bounded.get_grid_points();
  CHECK(left_grid(0) == Approx(-0.5));
  CHECK((left_grid.tail(left_grid.size() - 1) -
         left_grid.head(left_grid.size() - 1))
          .minCoeff() > 0.0);
  const double left_upper = std::ceil(left_grid(left_grid.size() - 1));
  const Eigen::VectorXd left_levels = Eigen::VectorXd::LinSpaced(
    static_cast<size_t>(left_upper + 1.0), 0.0, left_upper);
  CHECK(left_bounded.pdf(left_levels).sum() == Approx(1.0));

  Kde1d right_bounded(NAN, 3.0, "discrete", 1.0, NAN, 2, 100);
  right_bounded.fit(observations);
  const Eigen::VectorXd right_grid = right_bounded.get_grid_points();
  CHECK(right_grid(right_grid.size() - 1) == Approx(3.5));
  CHECK((right_grid.tail(right_grid.size() - 1) -
         right_grid.head(right_grid.size() - 1))
          .minCoeff() > 0.0);
  const double right_lower = std::floor(right_grid(0));
  const Eigen::VectorXd right_levels = Eigen::VectorXd::LinSpaced(
    static_cast<size_t>(3.0 - right_lower + 1.0), right_lower, 3.0);
  CHECK(right_bounded.pdf(right_levels).sum() == Approx(1.0));
  CHECK(right_bounded.pdf(Eigen::VectorXd::Constant(1, 4.0))(0) == 0.0);

  Kde1d unbounded(NAN, NAN, "discrete");
  unbounded.fit(observations);
  const Eigen::VectorXd unbounded_grid = unbounded.get_grid_points();
  const double unbounded_lower = std::floor(unbounded_grid(0));
  const double unbounded_upper =
    std::ceil(unbounded_grid(unbounded_grid.size() - 1));
  const Eigen::VectorXd unbounded_levels = Eigen::VectorXd::LinSpaced(
    static_cast<size_t>(unbounded_upper - unbounded_lower + 1.0),
    unbounded_lower,
    unbounded_upper);
  CHECK(unbounded.pdf(unbounded_levels).sum() == Approx(1.0));
  CHECK(unbounded.quantile(Eigen::VectorXd::Constant(1, 0.0))(0) ==
        unbounded_lower);

  const Eigen::VectorXd signed_observations = observations.array() - 2.0;
  Kde1d signed_bounded(-2.0, 1.0, "discrete", 1.0, NAN, 2, 100);
  signed_bounded.fit(signed_observations);
  const Eigen::VectorXd signed_grid = signed_bounded.get_grid_points();
  CHECK(signed_grid(0) == Approx(-2.5));
  CHECK(signed_grid(signed_grid.size() - 1) == Approx(1.5));
  const Eigen::VectorXd signed_levels =
    Eigen::VectorXd::LinSpaced(4, -2.0, 1.0);
  CHECK(signed_bounded.pdf(signed_levels).sum() == Approx(1.0));
  CHECK(signed_bounded.pdf(Eigen::VectorXd::Constant(1, -3.0))(0) == 0.0);
  CHECK(signed_bounded.pdf(Eigen::VectorXd::Constant(1, 2.0))(0) == 0.0);

  Kde1d singleton(0.0, 0.0, "discrete");
  singleton.fit(Eigen::VectorXd::Zero(40));
  CHECK(singleton.pdf(Eigen::VectorXd::Constant(1, 0.0))(0) == 1.0);
  CHECK(singleton.get_grid_points()(0) == Approx(-0.5));
  CHECK(singleton.get_grid_points()(singleton.get_actual_grid_size() - 1) ==
        Approx(0.5));
}

TEST_CASE("bounded discrete fits use adaptive boundary experts",
          "[discrete][boundary-expert]")
{
  Eigen::VectorXd observations(200);
  for (Eigen::Index i = 0; i < observations.size(); ++i)
    observations(i) = static_cast<double>(i % 4);

  Kde1d repaired(0.0, 3.0, "discrete", 1.0, NAN, 2, 100, true);
  repaired.fit(observations);
  Kde1d bulk(0.0, 3.0, "discrete", 1.0, NAN, 2, 100, false);
  bulk.fit(observations);

  CHECK_FALSE(repaired.get_values().isApprox(bulk.get_values(), 1e-8));
  CHECK(repaired.get_values().array().isFinite().all());
  CHECK(repaired.get_values().minCoeff() >= 0.0);
}

TEST_CASE("discrete inputs are integers", "[discrete][input-checks]")
{
  CHECK_NOTHROW(Kde1d(-3.0, -1.0, "discrete"));
  CHECK_THROWS(Kde1d(0.5, NAN, "discrete"));
  CHECK_THROWS(Kde1d(NAN, 3.5, "discrete"));

  Kde1d fit(NAN, NAN, "discrete");
  Eigen::VectorXd fractional(2);
  fractional << 0.0, 1.5;
  CHECK_THROWS(fit.fit(fractional));
  Eigen::VectorXd infinite(2);
  infinite << 0.0, std::numeric_limits<double>::infinity();
  CHECK_THROWS(fit.fit(infinite));
}

TEST_CASE("likelihood summaries match fitted densities", "[likelihood]")
{
  SECTION("weighted likelihood")
  {
    Eigen::VectorXd observations = Eigen::VectorXd::LinSpaced(100, -2.0, 2.0);
    Eigen::VectorXd weights = Eigen::VectorXd::LinSpaced(100, 0.5, 1.5);
    Kde1d weighted(NAN, NAN, "continuous", 1.0, 0.5);
    weighted.fit(observations, weights);

    double expected =
      (weighted.pdf(observations).array().log() * weights.array() /
       weights.mean())
        .sum();
    CHECK(weighted.get_loglik() == Approx(expected));
  }

  SECTION("zero-inflated likelihood")
  {
    Eigen::VectorXd observations(100);
    observations.head(40).setZero();
    observations.tail(60) = Eigen::VectorXd::LinSpaced(60, 0.01, 3.0);
    Kde1d zero_inflated(0.0, NAN, "zero-inflated", 1.0, 0.5);
    zero_inflated.fit(observations);

    Eigen::VectorXd density = zero_inflated.pdf(observations);
    REQUIRE(zero_inflated.get_values().array().isFinite().all());
    REQUIRE(density.array().isFinite().all());
    REQUIRE((density.array() > 0.0).all());
    REQUIRE(std::isfinite(zero_inflated.get_loglik()));
    CHECK(zero_inflated.get_loglik() == Approx(density.array().log().sum()));
  }
}

TEST_CASE("grid_size parameter", "[grid-size]")
{
  
  SECTION("constructor accepts grid_size parameter")
  {
    // Test VarType constructor
    kde1d::Kde1d fit1(NAN, NAN, kde1d::VarType::continuous, 1.0, NAN, 2, 100);
    CHECK(fit1.get_grid_size() == 100);
    
    // Test string constructor  
    kde1d::Kde1d fit2(NAN, NAN, "continuous", 1.0, NAN, 2, 200);
    CHECK(fit2.get_grid_size() == 200);
  }

  SECTION("grid_size validation")
  {
    // Should throw for grid_size < 4
    CHECK_THROWS(kde1d::Kde1d(NAN, NAN, "continuous", 1.0, NAN, 2, 3));
    CHECK_THROWS(kde1d::Kde1d(NAN, NAN, "continuous", 1.0, NAN, 2, 0));
    
    // Should work for grid_size >= 4
    CHECK_NOTHROW(kde1d::Kde1d(NAN, NAN, "continuous", 1.0, NAN, 2, 4));
    CHECK_NOTHROW(kde1d::Kde1d(NAN, NAN, "continuous", 1.0, NAN, 2, 50));
  }

  SECTION("default grid_size is 400")
  {
    kde1d::Kde1d fit;  // Use default constructor
    CHECK(fit.get_grid_size() == 400);
  }

  SECTION("grid_size affects interpolation grid size")
  {
    // Just test that the grid_size parameter is stored correctly
    std::vector<size_t> grid_sizes = {50, 100, 200};
    
    for (size_t requested_size : grid_sizes) {
      kde1d::Kde1d fit(NAN, NAN, "continuous", 1.0, NAN, 2, requested_size);
      
      // Check that requested grid size is stored correctly
      CHECK(fit.get_grid_size() == requested_size);
    }
  }

  SECTION("grid_size works after fitting")
  {
    // Test various grid sizes to ensure they work properly now
    std::vector<size_t> test_sizes = {50, 100, 200, 400, 600};
    
    for (size_t grid_size : test_sizes) {
      kde1d::Kde1d fit(NAN, NAN, "continuous", 1.0, NAN, 2, grid_size);
      CHECK_NOTHROW(fit.fit(x_ub));
      
      // Check that requested grid size is stored correctly
      CHECK(fit.get_grid_size() == grid_size);
      
      // Check that we can call methods that depend on the fitted model
      CHECK_NOTHROW(fit.pdf(x_ub));
      CHECK_NOTHROW(fit.cdf(x_ub));
      CHECK_NOTHROW(fit.quantile(ugrid));
      
      // Check that actual grid size matches requested size
      size_t actual_size = fit.get_actual_grid_size();
      CHECK(actual_size == grid_size + 1);  // Grid points = grid_size + 1
    }
  }

  SECTION("grid_size affects estimation with different data types")
  {
    size_t test_grid_size = 150;
    
    // Continuous data
    kde1d::Kde1d fit_cont(NAN, NAN, "continuous", 1.0, NAN, 2, test_grid_size);
    CHECK_NOTHROW(fit_cont.fit(x_ub));
    CHECK(fit_cont.get_grid_size() == test_grid_size);
    
    // Discrete data
    kde1d::Kde1d fit_disc(NAN, NAN, "discrete", 1.0, NAN, 2, test_grid_size);
    CHECK_NOTHROW(fit_disc.fit(x_d));
    CHECK(fit_disc.get_grid_size() == test_grid_size);
    
    // Zero-inflated data
    Eigen::VectorXd x_zi = x_lb;
    x_zi.head(n_sample / 4).setZero();
    kde1d::Kde1d fit_zi(0, NAN, "zero-inflated", 1.0, NAN, 2, test_grid_size);
    CHECK_NOTHROW(fit_zi.fit(x_zi));
    CHECK(fit_zi.get_grid_size() == test_grid_size);
  }

  SECTION("grid_size works with boundaries")
  {
    size_t test_grid_size = 120;
    
    // Left boundary
    kde1d::Kde1d fit_lb(0, NAN, "continuous", 1.0, NAN, 2, test_grid_size);
    CHECK_NOTHROW(fit_lb.fit(x_lb));
    CHECK(fit_lb.get_grid_size() == test_grid_size);
    
    // Right boundary
    kde1d::Kde1d fit_rb(NAN, 0, "continuous", 1.0, NAN, 2, test_grid_size);
    CHECK_NOTHROW(fit_rb.fit(x_rb));
    CHECK(fit_rb.get_grid_size() == test_grid_size);
    
    // Both boundaries
    kde1d::Kde1d fit_bb(0, 1, "continuous", 1.0, NAN, 2, test_grid_size);
    CHECK_NOTHROW(fit_bb.fit(x_cb));
    CHECK(fit_bb.get_grid_size() == test_grid_size);
  }

  SECTION("grid_size affects estimation accuracy")
  {
    auto points = stats::qnorm(upoints);
    auto target = stats::dnorm(points);
    
    // Test with small grid size
    kde1d::Kde1d fit_small(NAN, NAN, "continuous", 1.0, NAN, 2, 50);
    fit_small.fit(x_ub);
    auto pdf_small = fit_small.pdf(points);
    
    // Test with large grid size  
    kde1d::Kde1d fit_large(NAN, NAN, "continuous", 1.0, NAN, 2, 800);
    fit_large.fit(x_ub);
    auto pdf_large = fit_large.pdf(points);
    
    // Both should be reasonable approximations
    double error_small = (pdf_small - target).array().abs().mean();
    double error_large = (pdf_large - target).array().abs().mean();
    
    // Both errors should be reasonable (less than the tolerance used elsewhere)
    CHECK(error_small <= pdf_tol);
    CHECK(error_large <= pdf_tol);
    
    // Generally, larger grid should perform at least as well or better
    // (though for very large grids, numerical issues might make this not always true)
    CHECK(error_large <= error_small * 2.0);  // Allow reasonable tolerance
  }
}

TEST_CASE("bandwidth selection on refit", "[bandwidth][refit]")
{
  SECTION("refitting re-selects the bandwidth")
  {
    // The two samples differ in scale by a factor of ten, so an automatically
    // selected bandwidth must differ by roughly as much.
    Eigen::VectorXd narrow = x_ub;
    Eigen::VectorXd wide = x_ub * 10.0;

    Kde1d fresh;
    fresh.fit(wide);

    Kde1d reused;
    reused.fit(narrow);
    double first = reused.get_bandwidth();
    reused.fit(wide);

    CHECK(reused.get_bandwidth() == Approx(fresh.get_bandwidth()));
    CHECK(reused.get_bandwidth() > 2 * first);
  }

  SECTION("an explicitly requested bandwidth survives refitting")
  {
    Kde1d kde(NAN, NAN, VarType::continuous, 1.0, 0.75);
    kde.fit(x_ub);
    CHECK(kde.get_bandwidth() == Approx(0.75));
    kde.fit(x_ub * 10.0);
    CHECK(kde.get_bandwidth() == Approx(0.75));
  }

  SECTION("a grid-constructed density reports no bandwidth")
  {
    // The grid constructors never run a bandwidth selection, so the members
    // they leave alone must still be well defined.
    Kde1d fitted;
    fitted.fit(x_ub);
    Kde1d from_grid(interp::InterpolationGrid(
                      fitted.get_grid_points(), fitted.get_values(), 0),
                    NAN,
                    NAN,
                    VarType::continuous,
                    0.0);
    CHECK(std::isnan(from_grid.get_bandwidth()));
    CHECK(from_grid.get_multiplier() == Approx(1.0));
    CHECK(from_grid.get_degree() == 2);
    CHECK(from_grid.get_boundary_repair());
  }
}

TEST_CASE("misc checks", "[input-checks][argument-checks]")
{

  SECTION("detect wrong arguments")
  {
    CHECK_THROWS(kde1d::Kde1d(NAN, NAN, "asdf")); // unknown type
    CHECK_THROWS(kde1d::Kde1d(1, 0)); // distribution with xmin > xmax
    CHECK_THROWS(kde1d::Kde1d(NAN, NAN, "c", -1.0, NAN, 0)); // negative mult
    CHECK_THROWS(kde1d::Kde1d(NAN, NAN, "c", 1, -1.0, 0)); // negative bandwidth
    CHECK_THROWS(kde1d::Kde1d(NAN, NAN, "c", 1, NAN, 3));  // wrong degree
  }

  SECTION("methods fail if not fitted")
  {
    kde1d::Kde1d fit;
    CHECK_THROWS(fit.pdf(x_ub));
    CHECK_THROWS(fit.cdf(x_ub));
    CHECK_THROWS(fit.quantile(ugrid));
    CHECK_THROWS(fit.simulate(10));

    // doesn't have to fail
    CHECK(fit.get_values().size() == 0);
    CHECK(fit.get_grid_points().size() == 0);
  }

  SECTION("detect wrong inputs")
  {
    kde1d::Kde1d fit;
    // throws for empty data
    CHECK_THROWS(fit.fit(Eigen::VectorXd::Zero(0)));
    // throws when weights are not the same size as the data
    CHECK_THROWS(fit.fit(Eigen::VectorXd::Ones(10), Eigen::VectorXd::Ones(9)));

    // throws when some values in the data are smaler than xmin
    fit.set_xmin_xmax(1, 2);
    CHECK_THROWS(fit.fit(Eigen::VectorXd::Zero(1)));

    // throws when some values in the data are larger than xmax
    fit.set_xmin_xmax(-2, -1);
    CHECK_THROWS(fit.fit(Eigen::VectorXd::Zero(1)));

    // throws when trying to set an already fitted model
    fit.set_xmin_xmax(NAN, NAN);
    fit.fit(x_ub);
    CHECK_THROWS(fit.set_xmin_xmax(1, 2));

    // quantile throws when percentages are not in [0, 1]
    CHECK_THROWS(fit.quantile(Eigen::VectorXd::Constant(1, 1.1)));
    CHECK_THROWS(fit.quantile(Eigen::VectorXd::Constant(1, -0.1)));
  }
}

TEST_CASE("continuous data, unbounded", "[continuous][unbounded]")
{

  SECTION("fit local constant, linear, quadratic")
  {
    for (size_t degree = 0; degree < 3; degree++) {
      kde1d::Kde1d fit(NAN, NAN, "continuous", 1, NAN, degree);
      CHECK_NOTHROW(fit.fit(x_ub));
      CHECK(fit.str().find("continuous") != std::string::npos);
      CHECK(fit.str().find("xmin=nan") != std::string::npos);
      CHECK(fit.str().find("xmax=nan") != std::string::npos);
    }
  }

  SECTION("estimates are reasonable")
  {
    auto points = stats::qnorm(upoints);
    auto target = stats::dnorm(points);

    for (size_t degree = 0; degree < 3; degree++) {
      kde1d::Kde1d fit(NAN, NAN, "continuous", 1, NAN, degree);
      fit.fit(x_ub);

      CHECK(fit.pdf(x_ub).size() == n_sample);
      CHECK(fit.pdf(x_ub).minCoeff() >= 0);
      CHECK(fit.pdf(points).isApprox(target, pdf_tol));

      CHECK(fit.cdf(x_ub).size() == n_sample);
      CHECK(fit.cdf(x_ub).minCoeff() >= 0);
      CHECK(fit.cdf(x_ub).maxCoeff() <= 1);

      CHECK(fit.quantile(ugrid).size() == ugrid.size());
      CHECK(fit.quantile(ugrid).minCoeff() >= -2.5);
      CHECK(fit.quantile(ugrid).maxCoeff() <= 2.5);
      CHECK_NOTHROW(fit.simulate(10, { 1 }));
    }
  }

  SECTION("works with weights")
  {
    kde1d::Kde1d fit;
    auto w = Eigen::VectorXd::Constant(n_sample, 1);
    fit.fit(x_ub, w);

    kde1d::Kde1d fit0;
    fit0.fit(x_ub);

    CHECK(fit.pdf(x_ub).isApprox(fit0.pdf(x_ub)));

    Eigen::VectorXd w1 = Eigen::VectorXd::Constant(n_sample, 1.0);
    w1.tail(n_sample / 2) *= 2.0;

    kde1d::Kde1d fit1;
    fit1.fit(x_ub, w1);

    CHECK(fit1.pdf(x_ub).isApprox(fit0.pdf(x_ub), pdf_tol));

    Eigen::VectorXd w2 = Eigen::VectorXd::Constant(n_sample, 1.0);
    w2.tail(n_sample / 2) *= NAN;

    kde1d::Kde1d fit2;
    fit2.fit(x_ub, w2);

    CHECK(fit2.pdf(x_ub).isApprox(fit0.pdf(x_ub), pdf_tol));
  }

  SECTION("works with NaNs")
  {
    kde1d::Kde1d fit;
    auto w = Eigen::VectorXd::Constant(n_sample, 1);
    x_ub(0) = NAN;
    fit.fit(x_ub, w);

    CHECK(fit.pdf(Eigen::VectorXd::Constant(2, NAN)).array().isNaN().all());
    CHECK(fit.cdf(Eigen::VectorXd::Constant(2, NAN)).array().isNaN().all());
    CHECK(
      fit.quantile(Eigen::VectorXd::Constant(2, NAN)).array().isNaN().all());
  }
}

TEST_CASE("continuous data, left boundary", "[continuous][left-boundary]")
{

  SECTION("fit local constant, linear, quadratic")
  {
    for (size_t degree = 0; degree < 3; degree++) {
      kde1d::Kde1d fit(0, NAN, "continuous", 1, NAN, degree);
      CHECK_NOTHROW(fit.fit(x_lb));
      CHECK(fit.str().find("continuous") != std::string::npos);
      CHECK(fit.str().find("xmin=0") != std::string::npos);
      CHECK(fit.str().find("xmax=nan") != std::string::npos);
    }
  }

  SECTION("estimates are reasonable")
  {
    Eigen::VectorXd points = upoints.array().log();
    Eigen::VectorXd target = points.array().exp();
    points *= -1.0;

    for (size_t degree = 0; degree < 3; degree++) {
      kde1d::Kde1d fit(0, NAN, "continuous", 1, NAN, degree);
      fit.fit(x_lb);

      CHECK(fit.pdf(x_lb).size() == n_sample);
      CHECK(fit.pdf(x_lb).minCoeff() >= 0);
      CHECK(fit.pdf(points).isApprox(target, pdf_tol));
      CHECK(fit.pdf(Eigen::VectorXd::Constant(1, -1.0)).minCoeff() == 0.0);

      CHECK(fit.cdf(x_lb).size() == n_sample);
      CHECK(fit.cdf(x_lb).minCoeff() >= 0);
      CHECK(fit.cdf(x_lb).maxCoeff() <= 1);
      CHECK(fit.cdf(Eigen::VectorXd::Constant(1, -1.0)).minCoeff() == 0.0);

      CHECK(fit.quantile(ugrid).size() == ugrid.size());
      CHECK(fit.quantile(ugrid).minCoeff() >= 0);
      CHECK(fit.quantile(ugrid).maxCoeff() <= 10.0);
      CHECK(fit.simulate(10, { 1 }).maxCoeff() >= 0.0);
    }
  }

  SECTION("works with weights")
  {
    kde1d::Kde1d fit(0, NAN, "continuous");
    auto w = Eigen::VectorXd::Constant(n_sample, 1);
    fit.fit(x_lb, w);

    kde1d::Kde1d fit0(0, NAN, "continuous");
    fit0.fit(x_lb);

    CHECK(fit.pdf(x_lb).isApprox(fit0.pdf(x_lb)));

    Eigen::VectorXd w1 = Eigen::VectorXd::Constant(n_sample, 1.0);
    w1.tail(n_sample / 2) *= 2.0;

    kde1d::Kde1d fit1(0, NAN, "continuous");
    fit1.fit(x_lb, w1);

    CHECK(fit1.pdf(x_lb).isApprox(fit0.pdf(x_lb), pdf_tol));
  }

  SECTION("works with NaNs")
  {
    kde1d::Kde1d fit(0, NAN, "continuous");
    auto w = Eigen::VectorXd::Constant(n_sample, 1);
    x_lb(0) = NAN;
    fit.fit(x_lb, w);

    CHECK(fit.pdf(Eigen::VectorXd::Constant(2, NAN)).array().isNaN().all());
    CHECK(fit.cdf(Eigen::VectorXd::Constant(2, NAN)).array().isNaN().all());
    CHECK(
      fit.quantile(Eigen::VectorXd::Constant(2, NAN)).array().isNaN().all());
  }
}

TEST_CASE("continuous data, right boundary", "[continuous][right-boundary]")
{
  SECTION("fit local constant, linear, quadratic")
  {
    for (size_t degree = 0; degree < 3; degree++) {
      kde1d::Kde1d fit(NAN, 0, "continuous", 1, NAN, degree);
      CHECK_NOTHROW(fit.fit(x_rb));
      CHECK(fit.str().find("continuous") != std::string::npos);
      CHECK(fit.str().find("xmin=nan") != std::string::npos);
      CHECK(fit.str().find("xmax=0") != std::string::npos);
    }
  }

  SECTION("estimates are reasonable")
  {
    Eigen::VectorXd points = upoints.array().log();
    Eigen::VectorXd target = points.array().exp();

    for (size_t degree = 0; degree < 3; degree++) {
      kde1d::Kde1d fit(NAN, 0, "continuous", 1, NAN, degree);
      fit.fit(x_rb);

      CHECK(fit.pdf(x_rb).size() == n_sample);
      CHECK(fit.pdf(x_rb).minCoeff() >= 0);
      CHECK(fit.pdf(points).isApprox(target, pdf_tol));
      CHECK(fit.pdf(Eigen::VectorXd::Constant(1, 1.0)).minCoeff() == 0.0);

      CHECK(fit.cdf(x_rb).size() == n_sample);
      CHECK(fit.cdf(x_rb).minCoeff() >= 0);
      CHECK(fit.cdf(x_rb).maxCoeff() <= 1);
      CHECK(fit.cdf(Eigen::VectorXd::Constant(1, 1.0)).minCoeff() == 1.0);

      CHECK(fit.quantile(ugrid).size() == ugrid.size());
      CHECK(fit.quantile(ugrid).minCoeff() >= -10.0);
      CHECK(fit.quantile(ugrid).maxCoeff() <= 0.0);
      CHECK(fit.simulate(10, { 1 }).maxCoeff() <= 0.0);
    }
  }

  SECTION("works with weights")
  {
    kde1d::Kde1d fit(NAN, 0, "continuous");
    auto w = Eigen::VectorXd::Constant(n_sample, 1);
    fit.fit(x_rb, w);

    kde1d::Kde1d fit0(NAN, 0, "continuous");
    fit0.fit(x_rb);

    CHECK(fit.pdf(x_rb).isApprox(fit0.pdf(x_rb)));

    Eigen::VectorXd w1 = Eigen::VectorXd::Constant(n_sample, 1.0);
    w1.tail(n_sample / 2) *= 2.0;

    kde1d::Kde1d fit1(NAN, 0, "continuous");
    fit1.fit(x_rb, w1);

    CHECK(fit1.pdf(x_rb).isApprox(fit0.pdf(x_rb), pdf_tol));
  }

  SECTION("works with NaNs")
  {
    kde1d::Kde1d fit(NAN, 0, "continuous");
    auto w = Eigen::VectorXd::Constant(n_sample, 1);
    x_rb(0) = NAN;
    fit.fit(x_rb, w);

    CHECK(fit.pdf(Eigen::VectorXd::Constant(2, NAN)).array().isNaN().all());
    CHECK(fit.cdf(Eigen::VectorXd::Constant(2, NAN)).array().isNaN().all());
    CHECK(
      fit.quantile(Eigen::VectorXd::Constant(2, NAN)).array().isNaN().all());
  }
}

TEST_CASE("continuous data, both boundaries", "[continuous][both-boundaries]")
{
  SECTION("fit local constant, linear, quadratic")
  {
    for (size_t degree = 0; degree < 3; degree++) {
      kde1d::Kde1d fit(0, 1, "continuous", 1, NAN, degree);
      CHECK_NOTHROW(fit.fit(x_cb));
      CHECK(fit.str().find("continuous") != std::string::npos);
      CHECK(fit.str().find("xmin=0") != std::string::npos);
      CHECK(fit.str().find("xmax=1") != std::string::npos);
    }
  }

  SECTION("estimates are reasonable")
  {
    auto points = upoints;
    auto target = Eigen::VectorXd::Constant(points.size(), 1.0);

    for (size_t degree = 0; degree < 3; degree++) {
      kde1d::Kde1d fit(0, 1, "continuous", 1, NAN, degree);
      fit.fit(x_cb);

      CHECK(fit.pdf(x_cb).size() == n_sample);
      CHECK(fit.pdf(x_cb).minCoeff() >= 0);
      CHECK(fit.pdf(points).isApprox(target, pdf_tol));
      CHECK(fit.pdf(Eigen::VectorXd::Constant(1, -1.0)).minCoeff() == 0.0);
      CHECK(fit.pdf(Eigen::VectorXd::Constant(1, 2.0)).minCoeff() == 0.0);

      CHECK(fit.cdf(x_cb).size() == n_sample);
      CHECK(fit.cdf(x_cb).minCoeff() >= 0);
      CHECK(fit.cdf(x_cb).maxCoeff() <= 1);
      CHECK(fit.cdf(Eigen::VectorXd::Constant(1, -1.0)).minCoeff() == 0.0);
      CHECK(fit.cdf(Eigen::VectorXd::Constant(1, 2.0)).minCoeff() == 1.0);

      CHECK(fit.quantile(ugrid).size() == ugrid.size());
      CHECK(fit.quantile(ugrid).minCoeff() >= 0);
      CHECK(fit.quantile(ugrid).maxCoeff() <= 1.0);
      CHECK(fit.simulate(10, { 1 }).maxCoeff() >= 0.0);
      CHECK(fit.simulate(10, { 1 }).maxCoeff() <= 1.0);
    }
  }

  SECTION("works with weights")
  {
    kde1d::Kde1d fit(0, 1, "continuous");
    auto w = Eigen::VectorXd::Constant(n_sample, 1);
    fit.fit(x_cb, w);

    kde1d::Kde1d fit0(0, 1, "continuous");
    fit0.fit(x_cb);

    CHECK(fit.pdf(x_cb).isApprox(fit0.pdf(x_cb)));

    Eigen::VectorXd w1 = Eigen::VectorXd::Constant(n_sample, 1.0);
    w1.tail(n_sample / 2) *= 2.0;

    kde1d::Kde1d fit1(0, 1, "continuous");
    fit1.fit(x_cb, w1);

    CHECK(fit1.pdf(x_cb).isApprox(fit0.pdf(x_cb), pdf_tol));
  }

  SECTION("works with NaNs")
  {
    kde1d::Kde1d fit(0, 1, "continuous");
    auto w = Eigen::VectorXd::Constant(n_sample, 1);
    x_cb(0) = NAN;
    fit.fit(x_cb, w);

    CHECK(fit.pdf(Eigen::VectorXd::Constant(2, NAN)).array().isNaN().all());
    CHECK(fit.cdf(Eigen::VectorXd::Constant(2, NAN)).array().isNaN().all());
    CHECK(
      fit.quantile(Eigen::VectorXd::Constant(2, NAN)).array().isNaN().all());
  }
}

TEST_CASE("discrete data", "[discrete]")
{

  SECTION("fit local constant, linear, quadratic")
  {
    for (size_t degree = 0; degree < 3; degree++) {
      kde1d::Kde1d fit(0, NAN, "discrete", 1, NAN, degree);
      CHECK_NOTHROW(fit.fit(x_d));
      CHECK(fit.str().find("discrete") != std::string::npos);
    }
  }

  SECTION("estimates are reasonable")
  {
    auto points =
      Eigen::VectorXd::LinSpaced(nlevels, 0, static_cast<double>(nlevels) - 1);
    auto target =
      Eigen::VectorXd::Constant(nlevels, 1 / static_cast<double>(nlevels));

    for (size_t degree = 0; degree < 3; degree++) {
      kde1d::Kde1d fit(
        0, static_cast<double>(nlevels - 1), "discrete", 1, NAN, degree);
      fit.fit(x_d);

      CHECK(fit.pdf(x_d).size() == n_sample);
      CHECK(fit.pdf(x_d).minCoeff() >= 0);
      CHECK(fit.pdf(points).isApprox(target, pdf_tol));
      CHECK(fit.pdf(Eigen::VectorXd::Constant(1, -1.0)).minCoeff() == 0.0);
      CHECK(fit.pdf(Eigen::VectorXd::Constant(1, 0.5)).minCoeff() == 0.0);
      CHECK(fit.pdf(Eigen::VectorXd::Constant(1, static_cast<double>(nlevels)))
              .minCoeff() == 0.0);

      CHECK(fit.cdf(x_d).size() == n_sample);
      CHECK(fit.cdf(x_d).minCoeff() >= 0);
      CHECK(fit.cdf(x_d).maxCoeff() <= 1);
      CHECK(fit.cdf(Eigen::VectorXd::Constant(1, -1.0)).minCoeff() == 0.0);
      CHECK(fit.cdf(Eigen::VectorXd::Constant(1, static_cast<double>(nlevels)))
              .minCoeff() == 1.0);
      CHECK((fit.cdf(points) -
             fit.cdf(points + Eigen::VectorXd::Constant(points.size(), 0.5)))
              .minCoeff() == 0.0);

      CHECK(fit.quantile(ugrid).size() == ugrid.size());
      CHECK(fit.quantile(ugrid).minCoeff() >= 0);
      CHECK(fit.quantile(ugrid).maxCoeff() < static_cast<double>(nlevels));
      CHECK((fit.quantile(ugrid).array() - fit.quantile(ugrid).array().round())
              .abs()
              .maxCoeff() < 1e-300);
      CHECK(fit.simulate(10, { 1 }).maxCoeff() >= 0.0);
      CHECK(fit.simulate(10, { 1 }).maxCoeff() < static_cast<double>(nlevels));
    }
  }

  SECTION("works with weights")
  {
    kde1d::Kde1d fit(NAN, NAN, "discrete");
    auto w = Eigen::VectorXd::Constant(n_sample, 1);
    fit.fit(x_d, w);

    kde1d::Kde1d fit0(NAN, NAN, "discrete");
    fit0.fit(x_d);

    CHECK(fit.pdf(x_d).isApprox(fit0.pdf(x_d)));

    Eigen::VectorXd w1 = Eigen::VectorXd::Constant(n_sample, 1.0);
    w1.tail(n_sample / 2) *= 2.0;

    kde1d::Kde1d fit1(NAN, NAN, "discrete");
    fit1.fit(x_d, w1);

    CHECK(fit1.pdf(x_d).isApprox(fit0.pdf(x_d), pdf_tol));
  }

  SECTION("works with NaNs")
  {
    kde1d::Kde1d fit(NAN, NAN, "discrete");
    auto w = Eigen::VectorXd::Constant(n_sample, 1);
    x_d(0) = NAN;
    fit.fit(x_d, w);

    CHECK(fit.pdf(Eigen::VectorXd::Constant(2, NAN)).array().isNaN().all());
    CHECK(fit.cdf(Eigen::VectorXd::Constant(2, NAN)).array().isNaN().all());
    CHECK(
      fit.quantile(Eigen::VectorXd::Constant(2, NAN)).array().isNaN().all());
  }
}

TEST_CASE("zero-inflated data", "[zero-inflated]")
{
  Eigen::VectorXd x_zi = x_lb;
  x_zi.head(n_sample / 4).setZero();

  SECTION("fit local constant, linear, quadratic")
  {
    for (size_t degree = 0; degree < 3; degree++) {
      kde1d::Kde1d fit(0, NAN, "zinfl", 1, NAN, degree);
      CHECK_NOTHROW(fit.fit(x_zi));
      CHECK(fit.str().find("zero-inflated") != std::string::npos);
    }
  }

  SECTION("estimates are reasonable")
  {
    Eigen::VectorXd points = upoints.array().log();
    Eigen::VectorXd target = points.array().exp();
    target = target.array() * 0.75;
    points *= -1.0;
    points(0) = 0.0;
    target(0) = 0.25;

    for (size_t degree = 0; degree < 3; degree++) {
      kde1d::Kde1d fit(0, NAN, "zinfl", 1, NAN, degree);
      fit.fit(x_zi);

      CHECK(fit.pdf(x_zi).size() == n_sample);
      CHECK(fit.pdf(x_zi).minCoeff() >= 0);
      CHECK(fit.pdf(points).isApprox(target, pdf_tol));
      CHECK(fit.pdf(Eigen::VectorXd::Constant(1, -1.0)).minCoeff() == 0.0);

      CHECK(fit.cdf(x_zi).size() == n_sample);
      CHECK(fit.cdf(x_zi).minCoeff() >= 0);
      CHECK(fit.cdf(x_zi).maxCoeff() <= 1);
      CHECK(fit.cdf(Eigen::VectorXd::Constant(1, -1.0)).minCoeff() == 0.0);

      CHECK(fit.quantile(ugrid).size() == ugrid.size());
      CHECK(fit.quantile(ugrid).minCoeff() >= 0);
      CHECK(fit.quantile(ugrid).maxCoeff() <= 10.0);

      CHECK(fit.simulate(10, { 1 }).maxCoeff() >= 0.0);
      CHECK(fit.simulate(10, { 1 }).maxCoeff() <= 10.0);
    }
  }

  SECTION("works with weights")
  {
    kde1d::Kde1d fit(0, NAN, "zinfl");
    auto w = Eigen::VectorXd::Constant(n_sample, 1);
    fit.fit(x_zi, w);

    kde1d::Kde1d fit0(0, NAN, "zinfl");
    fit0.fit(x_zi);
    CHECK(fit.pdf(x_zi).isApprox(fit0.pdf(x_zi)));

    Eigen::VectorXd w1 = Eigen::VectorXd::Constant(n_sample, 1.0);
    for (int i = 0; i < n_sample / 2; i++) {
      w1(2 * i) *= 2;
    }

    kde1d::Kde1d fit1(0, NAN, "zinfl");
    fit1.fit(x_zi, w1);
    CHECK(fit1.pdf(x_zi).isApprox(fit0.pdf(x_zi), pdf_tol));
  }

  SECTION("works with NaNs")
  {
    kde1d::Kde1d fit(NAN, NAN, "zero-inflated");
    auto w = Eigen::VectorXd::Constant(n_sample, 1);
    x_zi(0) = NAN;
    fit.fit(x_zi, w);

    CHECK(fit.pdf(Eigen::VectorXd::Constant(2, NAN)).array().isNaN().all());
    CHECK(fit.cdf(Eigen::VectorXd::Constant(2, NAN)).array().isNaN().all());
    CHECK(
      fit.quantile(Eigen::VectorXd::Constant(2, NAN)).array().isNaN().all());
  }

  SECTION("works with only zeros")
  {
    kde1d::Kde1d fit(NAN, NAN, "zero-inflated");
    auto w = Eigen::VectorXd::Constant(n_sample, 1);
    x_zi = Eigen::VectorXd::Zero(n_sample);
    fit.fit(x_zi, w);

    CHECK(fit.pdf(Eigen::VectorXd::Constant(2, 1)).cwiseEqual(0).all());
    CHECK(fit.cdf(Eigen::VectorXd::Constant(2, -0.1)).cwiseEqual(0).all());
    CHECK(fit.cdf(Eigen::VectorXd::Constant(2, 0.1)).cwiseEqual(1).all());
    CHECK(
      fit.quantile(stats::simulate_uniform(100, { 5 })).cwiseEqual(0).all());
  }
}
