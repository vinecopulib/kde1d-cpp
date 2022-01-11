#include "../include/kde1d.hpp"

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

using namespace kde1d;

TEST_CASE("continuous data, unbounded", "[continuous][unbounded]")
{
  // continuous data
  Eigen::VectorXd x = stats::simulate_uniform(100, { 1 });

  SECTION("fit local constant")
  {
    kde1d::Kde1d fit(NAN, 1, NAN, NAN, 0);
    CHECK_NOTHROW(fit.fit(x));
  }

  SECTION("fit local linear")
  {
    kde1d::Kde1d fit(NAN, 1, NAN, NAN, 1);
    CHECK_NOTHROW(fit.fit(x));
  }

  SECTION("fit local quadratic")
  {
    kde1d::Kde1d fit;
    CHECK_NOTHROW(fit.fit(x));
  }

//   fit.pdf(x);
//   fit.cdf(x);
//   fit.quantile(x.cwiseMax(0));
}


TEST_CASE("discrete data", "[discrete]")
{
  // discrete data
  Eigen::VectorXi x = Eigen::VectorXi::LinSpaced(101, -50, 50);
  kde1d::Kde1d fit;
  fit.fit(x);
  // fit.pdf(x);
  // fit.cdf(x);
  // fit.quantile(Eigen::VectorXd::LinSpaced(100, 0.001, 0.999));
}
