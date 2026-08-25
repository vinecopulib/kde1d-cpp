#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch.hpp"
#include <kde1d.hpp>
#include <array>

TEST_CASE("continuous KDE performance", "[!benchmark]")
{
  Eigen::VectorXd unbounded = Eigen::VectorXd::LinSpaced(10000, -3.0, 3.0);
  unbounded.array() += 0.35 * (4.0 * unbounded.array()).sin();
  Eigen::VectorXd bounded =
    (1.0 / (1.0 + (-unbounded.array()).exp())).matrix();
  Eigen::VectorXd weights = Eigen::VectorXd::LinSpaced(10000, 0.5, 1.5);

  BENCHMARK("fit/end-to-end/bounded/tll/n=10000")
  {
    kde1d::Kde1d fit(0.0, 1.0, "continuous", 1.0, NAN, 2, 400);
    fit.fit(bounded);
    return fit.get_loglik();
  };

  BENCHMARK("fit/fixed-bandwidth/unbounded/tll/n=10000")
  {
    kde1d::Kde1d fit(NAN, NAN, "continuous", 1.0, 0.3, 2, 400);
    fit.fit(unbounded);
    return fit.get_loglik();
  };

  BENCHMARK("fit/fixed-bandwidth/bounded/tll/n=10000")
  {
    kde1d::Kde1d fit(0.0, 1.0, "continuous", 1.0, 0.05, 2, 400);
    fit.fit(bounded);
    return fit.get_loglik();
  };

  BENCHMARK("fit/fixed-bandwidth/bounded/tll/weighted/n=10000")
  {
    kde1d::Kde1d fit(0.0, 1.0, "continuous", 1.0, 0.05, 2, 400);
    fit.fit(bounded, weights);
    return fit.get_loglik();
  };

  kde1d::Kde1d unbounded_fit(
    NAN, NAN, "continuous", 1.0, 0.3, 2, 400);
  unbounded_fit.fit(unbounded);
  kde1d::Kde1d bounded_fit(0.0, 1.0, "continuous", 1.0, 0.05, 2, 400);
  bounded_fit.fit(bounded);
  Eigen::VectorXd bounded_grid_points = bounded_fit.get_grid_points();
  Eigen::VectorXd bounded_values = bounded_fit.get_values();
  kde1d::interp::InterpolationGrid bounded_grid(
    bounded_grid_points, bounded_values, 0);

  BENCHMARK("interpolation/construct/bounded/m=400")
  {
    return kde1d::interp::InterpolationGrid(
             bounded_grid_points, bounded_values, 0)
      .get_values()
      .sum();
  };

  for (Eigen::Index batch_size : std::array<Eigen::Index, 4>{ 10,
                                                              100,
                                                              1000,
                                                              10000 }) {
    Eigen::VectorXd unbounded_queries =
      Eigen::VectorXd::LinSpaced(batch_size, -4.0, 4.0);
    Eigen::VectorXd bounded_queries =
      Eigen::VectorXd::LinSpaced(batch_size, 1e-8, 1.0 - 1e-8);
    Eigen::VectorXd probabilities =
      Eigen::VectorXd::LinSpaced(batch_size, 1e-4, 1.0 - 1e-4);
    const std::string suffix = "/batch=" + std::to_string(batch_size);

    BENCHMARK("interpolation/interpolate/bounded" + suffix)
    {
      return bounded_grid.interpolate(bounded_queries).sum();
    };

    BENCHMARK("interpolation/construct-and-interpolate/bounded" + suffix)
    {
      kde1d::interp::InterpolationGrid grid(
        bounded_grid_points, bounded_values, 0);
      return grid.interpolate(bounded_queries).sum();
    };

    BENCHMARK("interpolation/integrate/bounded" + suffix)
    {
      return bounded_grid.integrate(bounded_queries, true).sum();
    };

    BENCHMARK("evaluate/pdf/unbounded/tll" + suffix)
    {
      return unbounded_fit.pdf(unbounded_queries).sum();
    };

    BENCHMARK("evaluate/pdf/bounded/tll" + suffix)
    {
      return bounded_fit.pdf(bounded_queries).sum();
    };

    BENCHMARK("evaluate/cdf/bounded/tll" + suffix)
    {
      return bounded_fit.cdf(bounded_queries).sum();
    };

    BENCHMARK("evaluate/quantile/bounded/tll" + suffix)
    {
      return bounded_fit.quantile(probabilities).sum();
    };
  }
}
