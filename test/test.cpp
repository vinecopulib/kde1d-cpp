#include "../include/kde1d.hpp"

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

using namespace kde1d;

Eigen::VectorXd ugrid = Eigen::VectorXd::LinSpaced(100, 0.001, 0.999);

TEST_CASE("continuous data, unbounded", "[continuous][unbounded]")
{
  // continuous data
  Eigen::VectorXd x = stats::simulate_uniform(100, { 1 });

  SECTION("fit local constant")
  {
    // no boundary
    kde1d::Kde1d fit(false, NAN, NAN, 1, NAN, 0);
    CHECK_NOTHROW(fit.fit(x));
  }

  SECTION("fit local linear")
  {
    kde1d::Kde1d fit(false, NAN, NAN, 1, NAN, 1);
    CHECK_NOTHROW(fit.fit(x));
  }

  SECTION("fit local quadratic")
  {
    kde1d::Kde1d fit;
    CHECK_NOTHROW(fit.fit(x));
  }

  SECTION("detect wrong arguments")
  {
    CHECK_THROWS(kde1d::Kde1d(true, 1, 0));
    CHECK_THROWS(kde1d::Kde1d(false, 1, 1, 1, NAN, 3));
  }

  SECTION("methods work")
  {
    kde1d::Kde1d fit;
    fit.fit(x);
    CHECK(fit.pdf(x).size() == x.size());
    CHECK(fit.pdf(x).minCoeff() >= 0);

    CHECK(fit.cdf(x).size() == x.size());
    CHECK(fit.cdf(x).minCoeff() >= 0);
    CHECK(fit.cdf(x).maxCoeff() <= 1);

    CHECK(fit.quantile(ugrid).size() == 100);
    CHECK(fit.quantile(ugrid).minCoeff() >= -1);
    CHECK(fit.quantile(ugrid).maxCoeff() <= 2);
  }

  SECTION("methods fail if not fitted")
  {
    kde1d::Kde1d fit;
    CHECK_THROWS(fit.pdf(x));
    CHECK_THROWS(fit.cdf(x));
    CHECK_THROWS(fit.quantile(ugrid));
    CHECK_THROWS(fit.simulate(10));

    // doesn't have to fail
    CHECK(fit.get_values().size() == 0);
    CHECK(fit.get_grid_points().size() == 0);
  }

  SECTION("estimates are reasonable")
  {
    x = stats::simulate_uniform(10000, { 1 });
    auto points = Eigen::VectorXd::LinSpaced(10, 0.1, 0.9);
    auto target = Eigen::VectorXd::Constant(10, 1.0);

    kde1d::Kde1d fit;
    fit.fit(x);
    CHECK(fit.pdf(points).isApprox(target, 0.2));

    fit = kde1d::Kde1d(false, NAN, NAN, 1, NAN, 0);
    fit.fit(x);
    CHECK(fit.pdf(points).isApprox(target, 0.2));

    fit = kde1d::Kde1d(false, NAN, NAN, 1, NAN, 1);
    fit.fit(x);
    CHECK(fit.pdf(points).isApprox(target, 0.2));
  }

  SECTION("works with weights")
  {
    kde1d::Kde1d fit;
    auto w = Eigen::VectorXd::Constant(x.size(), 1);
    fit.fit(x, w);

    kde1d::Kde1d fit0;
    fit0.fit(x);

    CHECK(fit.pdf(ugrid).isApprox(fit0.pdf(ugrid)));
  }
}

TEST_CASE("continuous data, left boundary", "[continuous][left-boundary]")
{
  // continuous data
  Eigen::VectorXd x = stats::simulate_uniform(100, { 1 });

  SECTION("fit local constant")
  {
    // no boundary
    kde1d::Kde1d fit(false, 0, NAN, 1, NAN, 0);
    CHECK_NOTHROW(fit.fit(x));
  }

  SECTION("fit local linear")
  {
    kde1d::Kde1d fit(false, 0, NAN, 1, NAN, 1);
    CHECK_NOTHROW(fit.fit(x));
  }

  SECTION("fit local quadratic")
  {
    kde1d::Kde1d fit(false, 0);
    CHECK_NOTHROW(fit.fit(x));
  }

  SECTION("methods work")
  {
    kde1d::Kde1d fit(false, 0);
    fit.fit(x);
    CHECK(fit.pdf(x).size() == x.size());
    CHECK(fit.pdf(x).minCoeff() >= 0);

    CHECK(fit.cdf(x).size() == x.size());
    CHECK(fit.cdf(x).minCoeff() >= 0);
    CHECK(fit.cdf(x).maxCoeff() <= 1);

    CHECK(fit.quantile(ugrid).size() == 100);
    CHECK(fit.quantile(ugrid).minCoeff() >= -1);
    CHECK(fit.quantile(ugrid).maxCoeff() <= 2);
  }

  SECTION("estimates are reasonable")
  {
    x = stats::simulate_uniform(10000, { 1 });
    auto points = Eigen::VectorXd::LinSpaced(10, 0.1, 0.9);
    auto target = Eigen::VectorXd::Constant(10, 1.0);

    kde1d::Kde1d fit(false, 0);
    fit.fit(x);
    CHECK(fit.pdf(points).isApprox(target, 0.2));

    fit = kde1d::Kde1d(false, 0, NAN, 1, NAN, 0);
    fit.fit(x);
    CHECK(fit.pdf(points).isApprox(target, 0.2));

    fit = kde1d::Kde1d(false, 0, NAN, 1, NAN, 1);
    fit.fit(x);
    CHECK(fit.pdf(points).isApprox(target, 0.2));
  }

  SECTION("works with weights")
  {
    kde1d::Kde1d fit(false, 0);
    auto w = Eigen::VectorXd::Constant(x.size(), 1);
    fit.fit(x, w);

    kde1d::Kde1d fit0(false, 0);
    fit0.fit(x);

    CHECK(fit.pdf(ugrid).isApprox(fit0.pdf(ugrid)));
  }
}

TEST_CASE("continuous data, right boundary", "[continuous][right-boundary]")
{
  // continuous data
  Eigen::VectorXd x = stats::simulate_uniform(100, { 1 });

  SECTION("fit local constant")
  {
    // no boundary
    kde1d::Kde1d fit(false, NAN, 1, 1, NAN, 0);
    CHECK_NOTHROW(fit.fit(x));
  }

  SECTION("fit local linear")
  {
    kde1d::Kde1d fit(false, NAN, 1, 1, NAN, 1);
    CHECK_NOTHROW(fit.fit(x));
  }

  SECTION("fit local quadratic")
  {
    kde1d::Kde1d fit(false, NAN, 1);
    CHECK_NOTHROW(fit.fit(x));
  }

  SECTION("methods work")
  {
    kde1d::Kde1d fit(false, NAN, 1);
    fit.fit(x);
    CHECK(fit.pdf(x).size() == x.size());
    CHECK(fit.pdf(x).minCoeff() >= 0);

    CHECK(fit.cdf(x).size() == x.size());
    CHECK(fit.cdf(x).minCoeff() >= 0);
    CHECK(fit.cdf(x).maxCoeff() <= 1);

    CHECK(fit.quantile(ugrid).size() == 100);
    CHECK(fit.quantile(ugrid).minCoeff() >= -1);
    CHECK(fit.quantile(ugrid).maxCoeff() <= 2);
  }

  SECTION("estimates are reasonable")
  {
    x = stats::simulate_uniform(10000, { 1 });
    auto points = Eigen::VectorXd::LinSpaced(10, 0.1, 0.9);
    auto target = Eigen::VectorXd::Constant(10, 1.0);

    kde1d::Kde1d fit(false, NAN, 1);
    fit.fit(x);
    CHECK(fit.pdf(points).isApprox(target, 0.2));

    fit = kde1d::Kde1d(false, NAN, 1, 1, NAN, 0);
    fit.fit(x);
    CHECK(fit.pdf(points).isApprox(target, 0.2));

    fit = kde1d::Kde1d(false, NAN, 1, 1, NAN, 1);
    fit.fit(x);
    CHECK(fit.pdf(points).isApprox(target, 0.2));
  }

  SECTION("works with weights")
  {
    kde1d::Kde1d fit(false, NAN, 1);
    auto w = Eigen::VectorXd::Constant(x.size(), 1);
    fit.fit(x, w);

    kde1d::Kde1d fit0(false, NAN, 1);
    fit0.fit(x);

    CHECK(fit.pdf(ugrid).isApprox(fit0.pdf(ugrid)));
  }
}

TEST_CASE("continuous data, both boundaries", "[continuous][both-boundaries]")
{
  // continuous data
  Eigen::VectorXd x = stats::simulate_uniform(100, { 1 });

  SECTION("fit local constant")
  {
    // no boundary
    kde1d::Kde1d fit(false, 0, 1, 1, NAN, 0);
    CHECK_NOTHROW(fit.fit(x));
  }

  SECTION("fit local linear")
  {
    kde1d::Kde1d fit(false, 0, 1, 1, NAN, 1);
    CHECK_NOTHROW(fit.fit(x));
  }

  SECTION("fit local quadratic")
  {
    kde1d::Kde1d fit(false, 0);
    CHECK_NOTHROW(fit.fit(x));
  }

  SECTION("methods work")
  {
    kde1d::Kde1d fit(false, 0, 1);
    fit.fit(x);
    CHECK(fit.pdf(x).size() == x.size());
    CHECK(fit.pdf(x).minCoeff() >= 0);

    CHECK(fit.cdf(x).size() == x.size());
    CHECK(fit.cdf(x).minCoeff() >= 0);
    CHECK(fit.cdf(x).maxCoeff() <= 1);

    CHECK(fit.quantile(ugrid).size() == 100);
    CHECK(fit.quantile(ugrid).minCoeff() >= -1);
    CHECK(fit.quantile(ugrid).maxCoeff() <= 2);
  }

  SECTION("estimates are reasonable")
  {
    x = stats::simulate_uniform(10000, { 1 });
    auto points = Eigen::VectorXd::LinSpaced(10, 0.1, 0.9);
    auto target = Eigen::VectorXd::Constant(10, 1.0);

    kde1d::Kde1d fit(false, 0, 1);
    fit.fit(x);
    CHECK(fit.pdf(points).isApprox(target, 0.2));

    fit = kde1d::Kde1d(false, 0, 1, 1, NAN, 0);
    fit.fit(x);
    CHECK(fit.pdf(points).isApprox(target, 0.2));

    fit = kde1d::Kde1d(false, 0, 1, 1, NAN, 1);
    fit.fit(x);
    CHECK(fit.pdf(points).isApprox(target, 0.2));
  }

  SECTION("works with weights")
  {
    kde1d::Kde1d fit(false, 0, 1);
    auto w = Eigen::VectorXd::Constant(x.size(), 1);
    fit.fit(x, w);

    kde1d::Kde1d fit0(false, 0, 1);
    fit0.fit(x);

    CHECK(fit.pdf(ugrid).isApprox(fit0.pdf(ugrid)));
  }
}

TEST_CASE("discrete data", "[discrete]")
{
  // discrete data
  Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(101, -50, 50);

  SECTION("fit local constant")
  {
    // no boundary
    kde1d::Kde1d fit(true, NAN, NAN, 1, NAN, 0);
    CHECK_NOTHROW(fit.fit(x));
  }

  SECTION("fit local linear")
  {
    kde1d::Kde1d fit(true, NAN, NAN, 1, NAN, 1);
    CHECK_NOTHROW(fit.fit(x));
  }

  SECTION("fit local quadratic")
  {
    kde1d::Kde1d fit(true);
    CHECK_NOTHROW(fit.fit(x));
  }

  SECTION("methods work")
  {
    kde1d::Kde1d fit(true);
    fit.fit(x);
    CHECK(fit.pdf(x).size() == x.size());
    CHECK(fit.pdf(x).minCoeff() >= 0);

    CHECK(fit.cdf(x).size() == x.size());
    CHECK(fit.cdf(x).minCoeff() >= 0);
    CHECK(fit.cdf(x).maxCoeff() <= 1);

    CHECK(fit.quantile(ugrid).size() == 100);
    CHECK(fit.quantile(ugrid).minCoeff() >= -50);
    CHECK(fit.quantile(ugrid).maxCoeff() <= 50);
    CHECK((fit.quantile(ugrid).array() - fit.quantile(ugrid).array().round())
            .abs()
            .maxCoeff() < 1e-300);
  }

  SECTION("estimates are reasonable")
  {
    x = stats::simulate_uniform(10000, { 1 });
    x = ((x.array() - 0.5) * 100).round();
    auto points = Eigen::VectorXd::LinSpaced(10, 0.1, 0.9);
    auto target = Eigen::VectorXd::Constant(10, 0.01);

    kde1d::Kde1d fit(true);
    fit.fit(x);
    CHECK(fit.pdf(points).isApprox(target, 0.2));

    fit = kde1d::Kde1d(true, NAN, NAN, 1, NAN, 0);
    fit.fit(x);
    CHECK(fit.pdf(points).isApprox(target, 0.2));

    fit = kde1d::Kde1d(true, NAN, NAN, 1, NAN, 1);
    fit.fit(x);
    CHECK(fit.pdf(points).isApprox(target, 0.2));
  }

  SECTION("works with weights")
  {
    kde1d::Kde1d fit(true);
    auto w = Eigen::VectorXd::Constant(x.size(), 1);
    fit.fit(x, w);

    kde1d::Kde1d fit0(true);
    fit0.fit(x);

    CHECK(fit.pdf(ugrid).isApprox(fit0.pdf(ugrid)));
  }
}
