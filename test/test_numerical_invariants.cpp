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
