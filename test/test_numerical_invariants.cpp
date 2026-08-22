#define CATCH_CONFIG_MAIN
#include "../include/kde1d.hpp"
#include "catch.hpp"

namespace {

void
check_continuous_invariants(const Eigen::VectorXd& observations,
                            double xmin,
                            double xmax,
                            size_t degree,
                            size_t grid_size,
                            bool weighted)
{
  INFO("degree=" << degree << ", grid_size=" << grid_size
                 << ", weighted=" << weighted);

  Eigen::VectorXd weights;
  if (weighted)
    weights = Eigen::VectorXd::LinSpaced(observations.size(), 0.5, 1.5);

  kde1d::Kde1d fit(xmin, xmax, "continuous", 1.0, 0.3, degree, grid_size);
  fit.fit(observations, weights);

  Eigen::VectorXd grid = fit.get_grid_points();
  Eigen::VectorXd density = fit.pdf(grid);
  Eigen::VectorXd distribution = fit.cdf(grid);

  REQUIRE(grid.size() >= static_cast<Eigen::Index>(grid_size));
  REQUIRE(fit.get_values().array().isFinite().all());
  REQUIRE(density.array().isFinite().all());
  REQUIRE(distribution.array().isFinite().all());
  CHECK(density.minCoeff() >= 0.0);
  CHECK(distribution.minCoeff() >= -1e-12);
  CHECK(distribution.maxCoeff() <= 1.0 + 1e-12);
  CHECK(distribution(0) == Approx(0.0).margin(1e-12));
  CHECK(distribution(distribution.size() - 1) == Approx(1.0).margin(1e-12));
  CHECK((distribution.tail(distribution.size() - 1) -
         distribution.head(distribution.size() - 1))
          .minCoeff() >= -1e-12);
  CHECK(std::isfinite(fit.get_loglik()));
  CHECK(std::isfinite(fit.get_edf()));

  if (!std::isnan(xmin)) {
    CHECK(fit.pdf(Eigen::VectorXd::Constant(1, xmin - 1e-8))(0) == 0.0);
    CHECK(fit.cdf(Eigen::VectorXd::Constant(1, xmin - 1e-8))(0) == 0.0);
  }
  if (!std::isnan(xmax)) {
    CHECK(fit.pdf(Eigen::VectorXd::Constant(1, xmax + 1e-8))(0) == 0.0);
    CHECK(fit.cdf(Eigen::VectorXd::Constant(1, xmax + 1e-8))(0) == 1.0);
  }
}

} // namespace

TEST_CASE("continuous fits satisfy numerical invariants",
          "[numerical-invariants]")
{
  Eigen::VectorXd unbounded = Eigen::VectorXd::LinSpaced(200, -2.0, 2.0);
  Eigen::VectorXd left_bounded = Eigen::VectorXd::LinSpaced(200, 0.01, 3.0);
  Eigen::VectorXd right_bounded = -left_bounded;
  Eigen::VectorXd bounded = Eigen::VectorXd::LinSpaced(200, 0.01, 0.99);

  for (size_t degree = 0; degree < 3; ++degree) {
    for (size_t grid_size : { 64, 400 }) {
      for (bool weighted : { false, true }) {
        check_continuous_invariants(
          unbounded, NAN, NAN, degree, grid_size, weighted);
        check_continuous_invariants(
          left_bounded, 0.0, NAN, degree, grid_size, weighted);
        check_continuous_invariants(
          right_bounded, NAN, 0.0, degree, grid_size, weighted);
        check_continuous_invariants(
          bounded, 0.0, 1.0, degree, grid_size, weighted);
      }
    }
  }
}

TEST_CASE("finite-support fits are scale equivariant",
          "[numerical-invariants][boundary-scale]")
{
  Eigen::VectorXd observations =
    Eigen::VectorXd::LinSpaced(2000, 0.5 / 2000.0, 1.0 - 0.5 / 2000.0);
  Eigen::VectorXd points(10);
  points << 0.0, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 0.1, 0.5, 0.9, 1.0;

  kde1d::Kde1d reference(0.0, 1.0, "continuous", 1.0, 0.3, 2, 400);
  reference.fit(observations);

  for (double scale : { 1e-4, 1e4 }) {
    INFO("scale=" << scale);
    kde1d::Kde1d scaled(
      0.0, scale, "continuous", 1.0, 0.3, 2, 400);
    scaled.fit(observations * scale);

    CHECK((scaled.get_grid_points() / scale)
            .isApprox(reference.get_grid_points(), 1e-11));
    CHECK((scaled.pdf(points * scale) * scale)
            .isApprox(reference.pdf(points), 1e-10));
    CHECK(scaled.cdf(points * scale).isApprox(reference.cdf(points), 1e-10));
  }
}

TEST_CASE("interpolation preserves nonuniform-grid reference values",
          "[interpolation][parity]")
{
  Eigen::VectorXd grid_points(7);
  grid_points << -2.0, -1.99, -1.2, -0.1, 0.0, 0.03, 2.5;
  Eigen::VectorXd values(7);
  values << 0.1, 0.8, 0.3, 1.2, 0.4, 1.1, 0.05;
  Eigen::VectorXd queries(14);
  queries << 2.5, -3.0, 0.015, -1.99, 1.1, -2.0, 3.0,
    -0.1, -1.5, 0.0, 0.03, -1.2, -1.2, NAN;

  kde1d::interp::InterpolationGrid grid(grid_points, values, 0);

  Eigen::VectorXd expected_interpolation(13);
  expected_interpolation << 0.05, 0.0, 0.72395374493927123, 0.8,
    8.6024753616471621, 0.1, 0.048985984426542499, 1.2,
    5.3485271720570431, 0.4, 1.1, 0.3, 0.3;
  Eigen::VectorXd interpolation = grid.interpolate(queries);
  CHECK(interpolation.head(13).isApprox(expected_interpolation, 1e-12));
  CHECK(std::isnan(interpolation(13)));

  Eigen::VectorXd expected_integral(13);
  expected_integral << 18.813458606991869, 0.0, 5.6615750998365577,
    0.0039240242616033786, 13.498277718530957, 0.0,
    18.813458606991869, 5.5893705472445951, 3.3631171509959734,
    5.6533162543153024, 5.6752953292140882, 4.0344600562080499,
    4.0344600562080499;
  Eigen::VectorXd integral = grid.integrate(queries);
  CHECK(integral.head(13).isApprox(expected_integral, 1e-12));
  CHECK(std::isnan(integral(13)));

  Eigen::VectorXd expected_normalized_integral(13);
  expected_normalized_integral << 1.0, 0.0, 0.30093217935656336,
    0.00020857537912487005, 0.71747986377764861, 0.0, 1.0,
    0.29709425916866511, 0.17876123796536264, 0.30049319332567026,
    0.30166145671402017, 0.21444542125330829, 0.21444542125330829;
  Eigen::VectorXd normalized_integral = grid.integrate(queries, true);
  CHECK(normalized_integral.head(13).isApprox(expected_normalized_integral,
                                               1e-12));
  CHECK(std::isnan(normalized_integral(13)));

  Eigen::VectorXd expected_normalized_values(7);
  expected_normalized_values << 0.005315343770062343,
    0.042522750160498744, 0.015946031310187028, 0.063784125240748113,
    0.021261375080249372, 0.058468781470685779,
    0.0026576718850311715;
  kde1d::interp::InterpolationGrid normalized(grid_points, values, 3);
  CHECK(normalized.get_values().isApprox(expected_normalized_values, 1e-12));

  grid.normalize(3);
  CHECK(grid.get_values().isApprox(expected_normalized_values, 1e-12));
  CHECK(grid.interpolate(queries).head(13).isApprox(
    expected_interpolation / expected_integral(0), 1e-12));
  CHECK(grid.integrate(queries).head(13).isApprox(
    expected_normalized_integral, 1e-12));
}
