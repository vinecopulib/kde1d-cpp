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
