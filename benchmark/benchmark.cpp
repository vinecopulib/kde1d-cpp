#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch.hpp"
#include <kde1d.hpp>

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
  Eigen::VectorXd unbounded_queries =
    Eigen::VectorXd::LinSpaced(100000, -4.0, 4.0);
  Eigen::VectorXd bounded_queries =
    Eigen::VectorXd::LinSpaced(100000, 1e-8, 1.0 - 1e-8);

  BENCHMARK("evaluate/pdf/unbounded/tll/n=100000")
  {
    return unbounded_fit.pdf(unbounded_queries).sum();
  };

  BENCHMARK("evaluate/pdf/bounded/tll/n=100000")
  {
    return bounded_fit.pdf(bounded_queries).sum();
  };

  BENCHMARK("evaluate/cdf/bounded/tll/n=100000")
  {
    return bounded_fit.cdf(bounded_queries).sum();
  };
}
