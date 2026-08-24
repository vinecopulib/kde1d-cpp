#pragma once

#include "dpik.hpp"
#include "interpolation.hpp"
#include "stats.hpp"
#include "tools.hpp"
#include <cmath>
#include <functional>

namespace kde1d {

enum class VarType
{
  continuous,
  discrete,
  zero_inflated
};

//! Local-polynomial density estimation in 1-d.
class Kde1d
{
public:
  // constructors
  Kde1d(double xmin,
        double xmax,
        VarType type,
        double multiplier = 1.0,
        double bandwidth = NAN,
        size_t degree = 2,
        size_t grid_size = 400);

  Kde1d(double xmin = NAN,
        double xmax = NAN,
        std::string type = "continuous",
        double multiplier = 1.0,
        double bandwidth = NAN,
        size_t degree = 2,
        size_t grid_size = 400);

  Kde1d(const interp::InterpolationGrid& grid,
        double xmin,
        double xmax,
        VarType type,
        double prob0_ = 0.0);

  Kde1d(const interp::InterpolationGrid& grid,
        double xmin = NAN,
        double xmax = NAN,
        std::string type = "continuous",
        double prob0_ = 0.0);

  void fit(const Eigen::VectorXd& x,
           const Eigen::VectorXd& weights = Eigen::VectorXd());

  // statistical functions
  Eigen::VectorXd pdf(const Eigen::VectorXd& x,
                      const bool& check_fitted = true) const;
  Eigen::VectorXd cdf(const Eigen::VectorXd& x,
                      const bool& check_fitted = true) const;
  Eigen::VectorXd quantile(const Eigen::VectorXd& x,
                           const bool& check_fitted = true) const;
  Eigen::VectorXd simulate(size_t n,
                           const std::vector<int>& seeds = {},
                           const bool& check_fitted = true) const;

  // getters
  //! @return the fitted density values on the interpolation grid.
  Eigen::VectorXd get_values() const { return grid_.get_values(); }
  //! @return the grid points used for interpolation (original scale).
  Eigen::VectorXd get_grid_points() const { return grid_.get_grid_points(); }
  //! @return the lower bound of the support (`NaN` if unbounded below).
  double get_xmin() const { return xmin_; }
  //! @return the upper bound of the support (`NaN` if unbounded above).
  double get_xmax() const { return xmax_; }
  //! @return the variable type (`continuous`, `discrete`, or `zero_inflated`).
  VarType get_type() const { return type_; }
  //! @return the variable type as the human-readable string
  //!   (`"continuous"`, `"discrete"`, or `"zero-inflated"`).
  std::string get_type_str() const { return this->as_str(type_); }
  //! @return the estimated point mass at zero (only used for the
  //!   `zero_inflated` type; `0` otherwise).
  double get_prob0() const { return prob0_; }
  //! @return the bandwidth multiplier supplied at construction.
  double get_multiplier() const { return multiplier_; }
  //! @return the bandwidth used by the most recent `fit()`, before the
  //!   multiplier is applied. Equals the bandwidth supplied at construction
  //!   when one was, and is `NaN` when none was and `fit()` has not run.
  double get_bandwidth() const { return bandwidth_; }
  //! @return the polynomial degree used by the local-likelihood
  //!   estimator (0, 1, or 2).
  size_t get_degree() const { return degree_; }
  //! @return the requested number of grid points (the value passed
  //!   to the constructor).
  size_t get_grid_size() const { return grid_size_; }
  //! @return the actual number of grid points after fitting (which
  //!   may differ slightly from `get_grid_size()` due to
  //!   boundary-snapping in `finalize_grid()`).
  size_t get_actual_grid_size() const { return grid_.get_grid_points().size(); }
  //! @return the effective degrees of freedom of the fitted estimator
  //!   (the sum of the per-observation influence values).
  double get_edf() const { return edf_; }
  //! @return the log-likelihood of the data under the fitted estimate.
  double get_loglik() const { return loglik_; }
  //! Updates the support bounds. Only valid before `fit()` has been called.
  //! @param xmin lower bound (`NaN` for unbounded).
  //! @param xmax upper bound (`NaN` for unbounded).
  void set_xmin_xmax(double xmin = NAN, double xmax = NAN);

  std::string str() const
  {
    std::stringstream ss;
    ss << "Kde1d("
       << "xmin=" << xmin_ << ", xmax=" << xmax_ << ", type='"
       << this->as_str(type_) << "'"
       << ", bandwidth=" << bandwidth_ << ", multiplier=" << multiplier_
       << ", degree=" << degree_ << ")";
    return ss.str();
  }

protected:
  void set_interpolation_grid(const interp::InterpolationGrid& grid);

private:
  struct BoundaryComponent
  {
    Eigen::VectorXd density;
    Eigen::VectorXd influence_num;
  };

  // data members
  interp::InterpolationGrid grid_;
  double xmin_;
  double xmax_;
  VarType type_;
  // Defaults matter: the grid constructors do not set these, so without them
  // the members are indeterminate for a density built from a grid.
  double multiplier_{ 1.0 };
  // The bandwidth as requested at construction (`NaN` for automatic
  // selection). Kept separate from `bandwidth_`, which holds the value the
  // most recent `fit()` settled on, so that refitting re-selects instead of
  // reusing the previous fit's bandwidth.
  double bandwidth_spec_{ NAN };
  double bandwidth_{ NAN };
  size_t degree_{ 2 };
  size_t grid_size_;
  double prob0_{ 0.0 };
  double loglik_{ NAN };
  double edf_{ NAN };
  // Regularizes one-sided power transforms in the units of the fitted data.
  double boundary_offset_{ NAN };
  // Makes one-sided transforms dimensionless and scale equivariant.
  double boundary_scale_{ 1.0 };
  static constexpr double K0_ = 0.3989425;

  // private methods
  void check_fitted() const;
  void check_notfitted() const;
  void check_xmin_xmax(const double& xmin, const double& xmax) const;
  void check_inputs(const Eigen::VectorXd& x,
                    const Eigen::VectorXd& weights = Eigen::VectorXd()) const;
  void check_boundaries(const Eigen::VectorXd& x) const;
  Eigen::VectorXd pdf_continuous(const Eigen::VectorXd& x) const;
  Eigen::VectorXd cdf_continuous(const Eigen::VectorXd& x) const;
  Eigen::VectorXd quantile_continuous(const Eigen::VectorXd& x) const;
  Eigen::VectorXd pdf_discrete(const Eigen::VectorXd& x) const;
  Eigen::VectorXd cdf_discrete(const Eigen::VectorXd& x) const;
  Eigen::VectorXd quantile_discrete(const Eigen::VectorXd& x) const;
  Eigen::VectorXd pdf_zi(const Eigen::VectorXd& x) const;
  Eigen::VectorXd cdf_zi(const Eigen::VectorXd& x) const;
  Eigen::VectorXd quantile_zi(const Eigen::VectorXd& x) const;

  Eigen::VectorXd kern_gauss(const Eigen::VectorXd& x);
  Eigen::MatrixXd fit_lp(const Eigen::VectorXd& x,
                         const Eigen::VectorXd& grid,
                         const Eigen::VectorXd& weights);
  double calculate_infl(const size_t& n,
                        const double& f0,
                        const double& f1,
                        const double& f2,
                        const double& bandwidth,
                        const double& s,
                        const double& weight);
  Eigen::VectorXd boundary_transform(const Eigen::VectorXd& x,
                                     bool inverse = false);
  Eigen::VectorXd boundary_correct(const Eigen::VectorXd& x,
                                   const Eigen::VectorXd& fhat);
  Eigen::VectorXd construct_grid_points(const Eigen::VectorXd& x);
  Eigen::VectorXd finalize_grid(Eigen::VectorXd& grid_points);
  double select_bandwidth(const Eigen::VectorXd& x,
                          double bandwidth,
                          double multiplier,
                          size_t degree,
                          const Eigen::VectorXd& weights) const;
  bool is_finite_endpoint(Eigen::VectorXd dist) const;
  BoundaryComponent fit_boundary_component(const Eigen::VectorXd& dist,
                                           const Eigen::VectorXd& eval_dist,
                                           double bandwidth) const;
  void repair_boundaries(const Eigen::VectorXd& x,
                         const Eigen::VectorXd& grid,
                         Eigen::VectorXd& influences);

  std::string as_str(VarType type) const;
  VarType as_enum(std::string type) const;
};

//! constructor for fitting the density estimate.
//! @param xmin lower bound for the support of the density, `NaN` means no
//!   boundary.
//! @param xmax upper bound for the support of the density, `NaN` means no
//!   boundary.
//! @param type variable type: `VarType::continuous`  for
//!   continuous variables, `VarType::discrete` for discrete integer
//!   variables, or `VarType::zero_inflated` for zero-inflated
//!   variables.
//! @param multiplier bandwidth multiplier (default is 1.0).
//! @param bandwidth positive bandwidth parameter (`NaN` means automatic
//! selection).
//! @param degree degree of the local polynomial.
//! @param grid_size number of grid points for the interpolation grid.
inline Kde1d::Kde1d(double xmin,
                    double xmax,
                    VarType type,
                    double multiplier,
                    double bandwidth,
                    size_t degree,
                    size_t grid_size)
  : xmin_(xmin)
  , xmax_(xmax)
  , type_(type)
  , multiplier_(multiplier)
  , bandwidth_spec_(bandwidth)
  , bandwidth_(bandwidth)
  , degree_(degree)
  , grid_size_(grid_size)
{
  this->check_xmin_xmax(xmin, xmax);
  if (multiplier <= 0.0) {
    throw std::invalid_argument("multiplier must be positive");
  }
  if (!std::isnan(bandwidth_) && (bandwidth_ <= 0.0)) {
    throw std::invalid_argument("bandwidth must be positive");
  }
  if (degree_ > 2) {
    throw std::invalid_argument("degree must be 0, 1 or 2");
  }
  if (grid_size_ < 4) {
    throw std::invalid_argument("grid_size must be at least 4");
  }
}

//! construct model from an already fit interpolation grid.
//! @param grid the interpolation grid.
//! @param xmin lower bound for the support of the density, `NaN` means no
//!   boundary.
//! @param xmax upper bound for the support of the density, `NaN` means no
//!   boundary.
//! @param type variable type: `VarType::continuous`  for
//!   continuous variables, `VarType::discrete` for discrete integer
//!   variables, or `VarType::zero_inflated` for zero-inflated
//!   variables.
//! @param prob0 point mass at 0.
inline Kde1d::Kde1d(const interp::InterpolationGrid& grid,
                    double xmin,
                    double xmax,
                    VarType type,
                    double prob0)
  : grid_(grid)
  , xmin_(xmin)
  , xmax_(xmax)
  , type_(type)
  , grid_size_(grid.get_grid_points().size())
  , prob0_(prob0)
{
  this->check_xmin_xmax(xmin, xmax);
  if ((prob0 < 0) || (prob0 > 1)) {
    throw std::invalid_argument("prob0 must lie in the interval [0, 1].");
  }
}

//! constructor for fitting the density estimate.
//! @param xmin lower bound for the support of the density, `NaN` means no
//!   boundary.
//! @param xmax upper bound for the support of the density, `NaN` means no
//!   boundary.
//! @param type variable type; must be one of {"c", "cont", "continuous"} for
//!   continuous variables, one of {"d", "disc", "discrete"} for discrete
//!   integer variables, or one of {"zi", "zinfl", "zero-inflated"} for
//!   zero-inflated variables.
//! @param multiplier bandwidth multiplier (default is 1.0).
//! @param bandwidth positive bandwidth parameter (`NaN` means automatic
//! selection).
//! @param degree degree of the local polynomial.
//! @param grid_size number of grid points for the interpolation grid.
inline Kde1d::Kde1d(double xmin,
                    double xmax,
                    std::string type,
                    double multiplier,
                    double bandwidth,
                    size_t degree,
                    size_t grid_size)
  : Kde1d(xmin,
          xmax,
          this->as_enum(type),
          multiplier,
          bandwidth,
          degree,
          grid_size)
{
}

//! construct model from an already fit interpolation grid.
//! @param grid the interpolation grid.
//! @param xmin lower bound for the support of the density, `NaN` means no
//!   boundary.
//! @param xmax upper bound for the support of the density, `NaN` means no
//!   boundary.
//! @param type variable type; must be one of {"c", "cont", "continuous"} for
//!   continuous variables, one of {"d", "disc", "discrete"} for discrete
//!   integer variables, or one of {"zi", "zinfl", "zero-inflated"} for
//!   zero-inflated variables.
//! @param prob0 point mass at 0.
inline Kde1d::Kde1d(const interp::InterpolationGrid& grid,
                    double xmin,
                    double xmax,
                    std::string type,
                    double prob0)
  : Kde1d(grid, xmin, xmax, this->as_enum(type), prob0)
{
}

//! Fits the kernel density estimate to data.
//! @param x vector of observations
//! @param weights vector of weights for each observation (optional).
inline void
Kde1d::fit(const Eigen::VectorXd& x, const Eigen::VectorXd& weights)
{
  check_inputs(x, weights);
  check_boundaries(x);

  // preprocessing for nans and jittering
  Eigen::VectorXd xx = x;
  Eigen::VectorXd w = weights;
  tools::remove_nans(xx, w);

  // Nonconstant case weights need separate expert semantics; constant weights
  // are equivalent to an unweighted sample.
  const bool unweighted = w.size() == 0 || w.minCoeff() == w.maxCoeff();
  const bool use_boundary_repair =
    type_ == VarType::continuous && xx.size() >= 16 && unweighted &&
    degree_ == 2 && (!std::isnan(xmin_) || !std::isnan(xmax_));
  Eigen::VectorXd boundary_observations;
  if (use_boundary_repair)
    boundary_observations = xx;

  if (w.size() > 0) {
    w /= w.mean();
  }

  if (type_ == VarType::zero_inflated) {
    if (w.size() == 0)
      w = Eigen::VectorXd::Ones(xx.size());
    w = (xx.array() == 0.0).select(Eigen::VectorXd::Zero(xx.size()), w);
    prob0_ = 1 - w.mean();
    xx =
      (w.array() == 0.0).select(Eigen::VectorXd::Constant(xx.size(), NAN), xx);
    tools::remove_nans(xx, w);
    if (xx.size() == 0) {
      bandwidth_ = NAN;
      loglik_ = 0.0;
      edf_ = 1.0;
      Eigen::VectorXd grid_points(5);
      grid_points << -2, -1, 0, 1, 2;
      auto values = Eigen::VectorXd::Constant(5, 0.0);
      grid_ = interp::InterpolationGrid(grid_points, values, 0);
      return;
    }
  } else if (type_ == VarType::discrete) {
    xx = stats::equi_jitter(xx);
  }

  if (type_ != VarType::discrete &&
      (std::isnan(xmin_) != std::isnan(xmax_))) {
    // Scaling the median boundary distance makes the transform equivariant
    // under changes of measurement units while remaining robust to outliers.
    Eigen::VectorXd distances;
    if (std::isnan(xmin_))
      distances = xmax_ - xx.array();
    else
      distances = xx.array() - xmin_;
    if (w.size() == 0) {
      boundary_scale_ = stats::median(distances);
    } else {
      boundary_scale_ = stats::quantile(
        distances, Eigen::VectorXd::Constant(1, 0.5), w)(0);
    }
    if (!(boundary_scale_ > 0.0))
      boundary_scale_ = distances.maxCoeff();
    if (!(boundary_scale_ > 0.0))
      boundary_scale_ = 1.0;
    boundary_offset_ = 1e-5 * boundary_scale_;
  }

  xx = boundary_transform(xx);

  // bandwidth selection (from the request, so refitting re-selects)
  bandwidth_ = select_bandwidth(xx, bandwidth_spec_, multiplier_, degree_, w);

  // fit model and evaluate in transformed domain
  Eigen::VectorXd grid_points = construct_grid_points(xx);
  Eigen::MatrixXd fitted = fit_lp(xx, boundary_transform(grid_points), w);

  // correct estimated density for transformation
  Eigen::VectorXd values = boundary_correct(grid_points, fitted.col(0));

  // order grid points from left to right
  grid_points = finalize_grid(grid_points);

  Eigen::VectorXd influences = fitted.col(1);
  if (std::isnan(xmin_) && !std::isnan(xmax_))
    influences.reverseInPlace();

  // construct interpolation grid
  // (3 iterations for normalization to a proper density)
  grid_ = interp::InterpolationGrid(grid_points, values, 3);
  if (use_boundary_repair)
    repair_boundaries(boundary_observations, grid_points, influences);

  // calculate log-likelihood of final estimate
  xx = boundary_transform(xx, true);
  if (type_ == VarType::discrete) {
    xx = xx.array().round();
  }
  
  if (w.size() == 0) {
    w = Eigen::VectorXd::Ones(xx.size());
  }

  loglik_ = (this->pdf(xx, false).array().log().array() * w.array()).sum();
  if (prob0_ > 0) {
    // For zero inflated data, all observations with value 0 have been removed,
    // so their likelihood contribution is missing. There were n * prob0_ such 
    // observations, each with log-likelihood contribution log(prob0_).
    loglik_ += static_cast<double>(x.size()) * prob0_ * std::log(prob0_);
  }

  // calculate effective degrees of freedom
  influences = influences.cwiseMin(3.0).cwiseMax(0);
  interp::InterpolationGrid infl_grid(grid_points, influences, 0);
  influences = infl_grid.interpolate(xx).array();
  edf_ = influences.sum() + static_cast<double>(prob0_ > 0);

  // store bandwidth in standardized format
  bandwidth_ = bandwidth_ / multiplier_;
}

//! Computes the pdf of the kernel density estimate by interpolation.
//! @param x vector of evaluation points.
//! @param check_fitted an optional logical to bypass the check.
//! @return a vector of pdf values.
inline Eigen::VectorXd
Kde1d::pdf(const Eigen::VectorXd& x, const bool& check_fitted) const
{
  if (check_fitted == true) {
    this->check_fitted();
  }
  check_inputs(x);

  switch (type_) {
    default:
      return pdf_continuous(x);
    case VarType::discrete:
      return pdf_discrete(x);
    case VarType::zero_inflated:
      return pdf_zi(x);
  }
}

inline Eigen::VectorXd
Kde1d::pdf_continuous(const Eigen::VectorXd& x) const
{
  Eigen::VectorXd fhat = grid_.interpolate(x);
  auto trunc = [](const double& xx) { return std::max(xx, 0.0); };
  fhat = tools::unaryExpr_or_nan(fhat, trunc);
  if (!std::isnan(xmin_))
    fhat = (x.array() < xmin_).select(0.0, fhat);
  if (!std::isnan(xmax_))
    fhat = (x.array() > xmax_).select(0.0, fhat);
  return fhat;
}

inline Eigen::VectorXd
Kde1d::pdf_discrete(const Eigen::VectorXd& x) const
{
  auto fhat = pdf_continuous(x);
  auto lb = std::floor(grid_.get_grid_min());
  auto ub = std::ceil(grid_.get_grid_max());
  Eigen::VectorXd lvs =
    Eigen::VectorXd::LinSpaced(static_cast<size_t>(ub - lb + 1), lb, ub);

  auto selected =
    (x.array() >= lb) && (x.array() <= ub) && (x.array() == x.array().round());
  fhat = fhat.array() * selected.cast<double>().array();

  // normalize
  fhat /= grid_.interpolate(lvs).sum();

  return fhat;
}

inline Eigen::VectorXd
Kde1d::pdf_zi(const Eigen::VectorXd& x) const
{
  auto ones = Eigen::VectorXd::Ones(x.size());
  return (x.array() == 0)
    .select(prob0_ * ones.array(), (1 - prob0_) * pdf_continuous(x).array());
}

//! Computes the cdf of the kernel density estimate by numerical
//! integration.
//! @param x vector of evaluation points.
//! @param check_fitted an optional logical to bypass the check.
//! @return a vector of cdf values.
inline Eigen::VectorXd
Kde1d::cdf(const Eigen::VectorXd& x, const bool& check_fitted) const
{
  if (check_fitted == true) {
    this->check_fitted();
  }
  check_inputs(x);

  switch (type_) {
    default:
      return cdf_continuous(x);
    case VarType::discrete:
      return cdf_discrete(x);
    case VarType::zero_inflated:
      return cdf_zi(x);
  }
}

inline Eigen::VectorXd
Kde1d::cdf_continuous(const Eigen::VectorXd& x) const
{
  return grid_.integrate(x, /* normalize */ true);
}

inline Eigen::VectorXd
Kde1d::cdf_discrete(const Eigen::VectorXd& x) const
{
  auto lb = std::floor(grid_.get_grid_min());
  auto ub = std::ceil(grid_.get_grid_max());
  Eigen::VectorXd lvs =
    Eigen::VectorXd::LinSpaced(static_cast<size_t>(ub - lb + 1), lb, ub);

  auto f_cum = pdf_discrete(lvs);
  for (Eigen::Index i = 1; i < f_cum.size(); ++i)
    f_cum(i) += f_cum(i - 1);

  return tools::unaryExpr_or_nan(x, [&](const double& xx) {
    if (xx < lb) {
      return 0.0;
    } else if (xx >= ub) {
      return 1.0;
    } else {
      return std::min(
        1.0, std::max(0.0, f_cum(static_cast<size_t>(xx - lb))));
    };
  });
}

inline Eigen::VectorXd
Kde1d::cdf_zi(const Eigen::VectorXd& x) const
{
  auto ones = Eigen::VectorXd::Ones(x.size());
  auto zeros = Eigen::VectorXd::Zero(x.size());
  Eigen::VectorXd zi = (x.array() >= 0).array().select(ones, zeros);
  return prob0_ * zi + (1 - prob0_) * (prob0_ < 1 ? cdf_continuous(x) : zeros);
}

//! Computes the quantile function by numerical inversion of the cdf.
//! @param x vector of evaluation points (probabilities in ``(0, 1)``).
//! @param check_fitted an optional logical to bypass the check.
//! @return a vector of quantiles.
inline Eigen::VectorXd
Kde1d::quantile(const Eigen::VectorXd& x, const bool& check_fitted) const
{
  if (check_fitted == true) {
    this->check_fitted();
  }
  if ((x.minCoeff() < 0) || (x.maxCoeff() > 1))
    throw std::invalid_argument("probabilities must lie in (0, 1).");

  switch (type_) {
    default:
      return quantile_continuous(x);
    case VarType::discrete:
      return quantile_discrete(x);
    case VarType::zero_inflated:
      return quantile_zi(x);
  }
}

inline Eigen::VectorXd
Kde1d::quantile_continuous(const Eigen::VectorXd& x) const
{
  return grid_.quantile(x);
}

inline Eigen::VectorXd
Kde1d::quantile_discrete(const Eigen::VectorXd& x) const
{
  auto lb = std::floor(grid_.get_grid_min());
  auto ub = std::ceil(grid_.get_grid_max());
  auto nlevels = static_cast<size_t>(ub - lb + 1);
  Eigen::VectorXd lvs = Eigen::VectorXd::LinSpaced(nlevels, lb, ub);

  auto p = cdf_discrete(lvs);
  auto quan = [&](const double& pp) {
    const double* level = std::lower_bound(p.data(), p.data() + p.size(), pp);
    return lvs(std::min(static_cast<Eigen::Index>(level - p.data()),
                        lvs.size() - 1));
  };

  return tools::unaryExpr_or_nan(x, quan);
}

inline Eigen::VectorXd
Kde1d::quantile_zi(const Eigen::VectorXd& x) const
{
  if (prob0_ >= 1.0)
    return tools::unaryExpr_or_nan(x, [](const double&) { return 0.0; });

  Eigen::VectorXd qs(x.size());
  auto p0 = this->cdf(Eigen::VectorXd::Zero(1), false)(0);
  auto newx = (x.array() <= p0 - prob0_)
                .select(x / (1 - prob0_),
                        (x.array() - prob0_).cwiseMax(0.0) / (1 - prob0_));
  qs = this->quantile_continuous(newx);
  for (Eigen::Index i = 0; i < x.size(); i++) {
    if ((x(i) > p0 - prob0_) && (x(i) <= p0)) {
      qs(i) = 0;
    }
  }
  return qs;
}

//! Simulates data from the fitted density.
//! @param n the number of observations to simulate.
//! @param seeds an optional vector of seeds.
//! @param check_fitted an optional logical to bypass the check.
//! @return simulated observations from the kernel density.
inline Eigen::VectorXd
Kde1d::simulate(size_t n,
                const std::vector<int>& seeds,
                const bool& check_fitted) const
{
  if (check_fitted == true) {
    this->check_fitted();
  }
  auto u = stats::simulate_uniform(n, seeds);
  return this->quantile(u);
}

//! Gaussian kernel (truncated at +/- 5).
//! @param x vector of evaluation points.
inline Eigen::VectorXd
Kde1d::kern_gauss(const Eigen::VectorXd& x)
{
  auto f = [](double xx) {
    // truncate at +/- 5
    if (std::fabs(xx) > 5.0)
      return 0.0;
    // otherwise calculate normal pdf (orrect for truncation)
    return stats::dnorm(Eigen::VectorXd::Constant(1, xx))(0) / 0.999999426;
  };
  return x.unaryExpr(f);
}

//! (analytically) evaluates the kernel density estimate and its influence
//! function on a user-supplied grid.
//! @param x_ev evaluation points.
//! @param x observations.
//! @param weights vector of weights for each observation (can be empty).
//! @return a two-column matrix containing the density estimate in the first
//!   and the influence function in the second column.
inline Eigen::MatrixXd
Kde1d::fit_lp(const Eigen::VectorXd& x,
              const Eigen::VectorXd& grid_points,
              const Eigen::VectorXd& weights)
{
  size_t m = grid_points.size();
  fft::KdeFFT kde_fft(
    x, bandwidth_, grid_points(0), grid_points(m - 1), weights, m - 1);
  Eigen::VectorXd f0 = kde_fft.kde_drv(0);
  Eigen::VectorXd f1(f0.size()), f2(f0.size());
  const double density_floor =
    std::numeric_limits<double>::epsilon() * f0.cwiseAbs().maxCoeff();

  Eigen::VectorXd wbin = Eigen::VectorXd::Ones(m);
  if (weights.size()) {
    // compute the average weight per cell
    auto wcount = kde_fft.get_bin_counts();
    auto count = tools::linbin(x,
                               grid_points(0),
                               grid_points(m - 1),
                               m - 1,
                               Eigen::VectorXd::Ones(x.size()));
    wbin = (count.array() == 0).select(
      Eigen::VectorXd::Zero(count.size()),
      wcount.cwiseQuotient(count)
    );
  }

  Eigen::MatrixXd res(f0.size(), 2);
  res.col(0) = f0;
  res.col(1) =
    K0_ / (static_cast<double>(x.size()) * bandwidth_) * wbin.cwiseQuotient(f0);
  if (degree_ == 0)
    return res;

  // degree > 0
  f1 = kde_fft.kde_drv(1);
  Eigen::VectorXd S = Eigen::VectorXd::Constant(f0.size(), bandwidth_);
  Eigen::VectorXd b = f1.cwiseQuotient(f0);
  if (degree_ == 2) {
    f2 = kde_fft.kde_drv(2);
    // D/R is notation from Hjort and Jones' AoS paper
    Eigen::VectorXd D = f2.cwiseQuotient(f0) - b.cwiseProduct(b);
    Eigen::VectorXd R = 1 / (1.0 + bandwidth_ * bandwidth_ * D.array()).sqrt();
    // this is our notation
    S = (R / bandwidth_).array().pow(2);
    b *= bandwidth_ * bandwidth_;
    res.col(0) = bandwidth_ * S.cwiseSqrt().cwiseProduct(res.col(0));
  }
  res.col(0) = res.col(0).array() * (-0.5 * b.array().pow(2) * S.array()).exp();

  for (size_t k = 0; k < m; k++) {
    if (!std::isfinite(f0(k)) || f0(k) <= density_floor) {
      res.row(k).setZero();
      continue;
    }
    res(k, 1) =
      calculate_infl(x.size(), f0(k), f1(k), f2(k), bandwidth_, S(k), wbin(k));
    if (!std::isfinite(res(k, 0)) || res(k, 0) < 0.0)
      res.row(k).setZero();
  }

  return res;
}

//! calculate influence for data point for density estimate based on
//! quantities pre-computed in `fit_lp()`.
inline double
Kde1d::calculate_infl(const size_t& n,
                      const double& f0,
                      const double& f1,
                      const double& f2,
                      const double& bandwidth,
                      const double& s,
                      const double& weight)
{
  double M_inverse00;
  double B = bandwidth * bandwidth;
  if (degree_ == 0) {
    M_inverse00 = 1 / f0;
  } else if (degree_ == 1) {
    Eigen::Matrix2d M;
    M(0, 0) = f0;
    M(0, 1) = B * f1;
    M(1, 0) = M(0, 1);
    M(1, 1) = f0 * B + B * f1 * f1 * B / f0;
    M_inverse00 = M.inverse()(0, 0);
  } else {
    Eigen::Matrix3d M;
    M(0, 0) = f0;
    M(0, 1) = B * f1;
    M(1, 0) = M(0, 1);
    M(1, 1) = B * f2 * B + B * f0;
    M(2, 0) = M(1, 1) / 2;
    M(0, 2) = M(1, 1) / 2;
    double s2 = B * f1 / f0;
    M(1, 2) = f0 / 2 * (3 / s * s2 + std::pow(s2, 3));
    M(2, 1) = M(1, 2);
    M(2, 2) =
      f0 / 4 * (3 / (s * s) + 6 / s * std::pow(s2, 2) + std::pow(s2, 4));
    M_inverse00 = M.inverse()(0, 0);
  }

  return K0_ * weight / (static_cast<double>(n) * bandwidth) * M_inverse00;
}

//! transformations for density estimates with bounded support.
//! @param x evaluation points.
//! @param inverse whether the inverse transformation should be applied.
//! @return the transformed evaluation points.
inline Eigen::VectorXd
Kde1d::boundary_transform(const Eigen::VectorXd& x, bool inverse)
{
  if (type_ == VarType::discrete) {
    return x; // no transform for discrete variables
  }

  Eigen::VectorXd x_new = x;
  if (!inverse) {
    if (!std::isnan(xmin_) && !std::isnan(xmax_)) {
      // two boundaries -> probit transform
      auto rng = xmax_ - xmin_;
      x_new = (x.array() - xmin_ + 5e-5 * rng) / (1.0001 * rng);
      x_new = stats::qnorm(x_new);
    } else if (!std::isnan(xmin_)) {
      // left boundary -> power-3/4 transform
      x_new = 4.0 * (((boundary_offset_ + x.array() - xmin_) /
                      boundary_scale_)
                       .pow(0.25) -
                     std::pow(boundary_offset_ / boundary_scale_, 0.25));
    } else if (!std::isnan(xmax_)) {
      // right boundary -> reflected power-3/4 transform
      x_new = 4.0 * (((boundary_offset_ + xmax_ - x.array()) /
                      boundary_scale_)
                       .pow(0.25) -
                     std::pow(boundary_offset_ / boundary_scale_, 0.25));
    } else {
      // no boundary -> no transform
    }
  } else {
    if (!std::isnan(xmin_) && !std::isnan(xmax_)) {
      // two boundaries -> probit transform
      auto rng = xmax_ - xmin_;
      x_new = stats::pnorm(x).array() * 1.0001 * rng + xmin_ - 5e-5 * rng;
    } else if (!std::isnan(xmin_)) {
      // left boundary -> inverse power-3/4 transform
      x_new = boundary_scale_ *
                (x.array() / 4.0 +
                 std::pow(boundary_offset_ / boundary_scale_, 0.25))
                  .pow(4) +
              xmin_ - boundary_offset_;
    } else if (!std::isnan(xmax_)) {
      // right boundary -> inverse reflected power-3/4 transform
      x_new = xmax_ + boundary_offset_ -
              boundary_scale_ *
                (x.array() / 4.0 +
                 std::pow(boundary_offset_ / boundary_scale_, 0.25))
                  .pow(4);
    } else {
      // no boundary -> no transform
    }
  }

  return x_new;
}

//! corrects the density estimate for a preceding boundary transformation of
//! the data.
//! @param x evaluation points (in original domain).
//! @param fhat the density estimate evaluated in the transformed domain.
//! @return corrected density estimates at `x`.
inline Eigen::VectorXd
Kde1d::boundary_correct(const Eigen::VectorXd& x, const Eigen::VectorXd& fhat)
{
  if (type_ == VarType::discrete) {
    return fhat; // no transform for discrete variables
  }

  Eigen::VectorXd corr_term(fhat.size());
  if (!std::isnan(xmin_) && !std::isnan(xmax_)) {
    // two boundaries -> probit transform
    auto rng = xmax_ - xmin_;
    corr_term = (x.array() - xmin_ + 5e-5 * rng) / (xmax_ - xmin_ + 1e-4 * rng);
    corr_term = stats::dnorm(stats::qnorm(corr_term));
    corr_term /= (xmax_ - xmin_ + 1e-4 * rng);
    corr_term = 1.0 / corr_term.array();
  } else if (!std::isnan(xmin_)) {
    // left boundary -> power-3/4 transform
    corr_term = ((boundary_offset_ + x.array() - xmin_) / boundary_scale_)
                  .pow(-0.75) /
                boundary_scale_;
  } else if (!std::isnan(xmax_)) {
    // right boundary -> reflected power-3/4 transform
    corr_term = ((boundary_offset_ + xmax_ - x.array()) / boundary_scale_)
                  .pow(-0.75) /
                boundary_scale_;
  } else {
    // no boundary -> no transform
    corr_term.fill(1.0);
  }

  Eigen::VectorXd f_corr = fhat.cwiseProduct(corr_term);
  if (std::isnan(xmin_) && !std::isnan(xmax_))
    f_corr.reverseInPlace();

  return f_corr;
}

//! constructs a grid later used for interpolation
//! @param x vector of observations.
//! @return a grid of size 400.
inline Eigen::VectorXd
Kde1d::construct_grid_points(const Eigen::VectorXd& x)
{
  Eigen::VectorXd rng(2);
  rng << x.minCoeff(), x.maxCoeff();
  if (type_ == VarType::discrete) {
    // Discrete estimates use jittered observations without transformation.
  } else if (std::isnan(xmin_) && std::isnan(xmax_)) {
    rng(0) -= 4 * bandwidth_;
    rng(1) += 4 * bandwidth_;
  } else if (!std::isnan(xmin_) && !std::isnan(xmax_)) {
    Eigen::VectorXd boundaries(2);
    boundaries << xmin_, xmax_;
    rng = boundary_transform(boundaries);
  } else {
    rng(0) = boundary_transform(Eigen::VectorXd::Constant(
      1, std::isnan(xmin_) ? xmax_ : xmin_))(0);
    rng(1) += 4 * bandwidth_;
  }
  auto zgrid = Eigen::VectorXd::LinSpaced(grid_size_ + 1, rng(0), rng(1));
  Eigen::VectorXd grid_points = boundary_transform(zgrid, true);

  // Avoid round-off at finite support boundaries before evaluating the fit.
  if (type_ != VarType::discrete) {
    if (!std::isnan(xmin_))
      grid_points(0) = xmin_;
    if (!std::isnan(xmax_))
      grid_points(std::isnan(xmin_) ? 0 : grid_size_) = xmax_;
  }

  return grid_points;
}

//! orders grid points from left to right.
//! @param grid_points the grid points.
inline Eigen::VectorXd
Kde1d::finalize_grid(Eigen::VectorXd& grid_points)
{
  if (std::isnan(xmin_) && !std::isnan(xmax_))
    grid_points.reverseInPlace();
  if (type_ == VarType::discrete) {
    if (!std::isnan(xmin_))
      grid_points(0) = xmin_;
    if (!std::isnan(xmax_))
      grid_points(grid_points.size() - 1) = xmax_;
  }

  return grid_points;
}

//! Classifies a support endpoint from its ordered distances @f$d@f$ to the
//! endpoint. For @f$k = \min(n-1, \lceil 2\sqrt n \rceil)@f$, the lower-tail
//! index is
//! @f[
//! \widehat\beta = \frac{k}{\sum_{i=1}^k
//!   \log\{d_{(k+1)} / d_{(i)}\}}.
//! @f]
//! It is treated as finite when @f$\widehat\beta \geq 0.9@f$ and its
//! one-sided 95% lower confidence bound does not exceed one. Ambiguous and
//! numerically degenerate cases retain the bulk fit.
//! @param dist distances of the observations from the support endpoint; copied
//!   because the classifier sorts them.
//! @return whether the endpoint is confidently classified as finite.
inline bool
Kde1d::is_finite_endpoint(Eigen::VectorXd dist) const
{
  std::sort(dist.data(), dist.data() + dist.size());

  // k is the number of lower order statistics; dist_k1 is d_(k+1).
  const size_t k =
    std::min(static_cast<size_t>(dist.size() - 1),
             static_cast<size_t>(std::ceil(2.0 * std::sqrt(dist.size()))));
  const double dist_k1 = dist(k);
  if (!(dist_k1 > 0.0) || !std::isfinite(dist_k1))
    return false;

  // denom is the denominator of the lower-tail index beta = k / denom.
  double denom = 0.0;
  const double dist_min = std::numeric_limits<double>::epsilon() * dist_k1;
  for (size_t i = 0; i < k; ++i) {
    denom += std::log(dist_k1 / std::max(dist(i), dist_min));
  }
  if (!(denom > 0.0) || !std::isfinite(denom))
    return false;

  const double beta = static_cast<double>(k) / denom;
  // 1.64485 is the 95% standard-normal quantile used by the R selector.
  return beta >= 0.9 && beta * (1.0 - 1.6448536269514722 / std::sqrt(k)) <= 1.0;
}

//! Fits the average of Gaussian local-linear and local-quadratic equivalent
//! kernels from sample distances @f$d@f$ to endpoint-distance grid @f$t@f$.
//! For degree @f$p@f$, shared bandwidth @f$h@f$, truncated Gaussian moments
//! @f$\mu_r(a)@f$, and @f$M_p(a)=[\mu_{j+k}(a)]_{j,k=0}^p@f$,
//! @f[
//! \begin{align*}
//! a &= t/h, & u_i &= (t-d_i)/h, \\[2pt]
//! \widehat f_p(t;h)
//!   &= \frac{1}{nh}\sum_{i=1}^n \phi(u_i)\,
//!      e_0^\mathsf{T}M_p(a)^{-1}
//!      (1,u_i,\ldots,u_i^p)^\mathsf{T}, \\[2pt]
//! \widehat f_\partial(t)
//!   &= \{\widehat f_1(t;h)+\widehat f_2(t;h)\}/2.
//! \end{align*}
//! @f]
//! Writing @f$f^{(r)}@f$ for the @f$r@f$th derivative of the ordinary
//! Gaussian KDE and @f$c_{pj}@f$ for the entries of the first row of
//! @f$M_p(a)^{-1}@f$, the same estimator is
//! @f[
//! \begin{align*}
//! \widehat f_\partial(t)
//!   &= \frac{1}{2}\{(c_{10}+c_{20}+c_{22})f(t)
//!     -h(c_{11}+c_{21})f^{(1)}(t)+h^2c_{22}f^{(2)}(t)\}.
//! \end{align*}
//! @f]
//! This identity evaluates the three convolutions by FFT.
//! Observations contribute equally to both kernels. The returned influence
//! numerator is the corresponding average diagonal kernel contribution.
//! @param dist ordered distances of the observations from the support endpoint.
//! @param eval_dist endpoint distances at which to evaluate the component.
//! @param bandwidth shared bandwidth of both equivalent kernels.
//! @return the density and influence numerator on `eval_dist`.
inline Kde1d::BoundaryComponent
Kde1d::fit_boundary_component(const Eigen::VectorXd& dist,
                              const Eigen::VectorXd& eval_dist,
                              double bandwidth) const
{
  const double h = bandwidth;
  if (!(h > 0.0) || !std::isfinite(h))
    return {};

  BoundaryComponent fit{ Eigen::VectorXd::Zero(eval_dist.size()),
                         Eigen::VectorXd::Zero(eval_dist.size()) };
  if (eval_dist.size() == 0)
    return fit;

  // Compute the ordinary Gaussian KDE and its first two derivatives on a
  // regular grid, then interpolate them to the endpoint grid.
  const double fft_upper = eval_dist(eval_dist.size() - 1) + 6.0 * h;
  const double* last =
    std::upper_bound(dist.data(), dist.data() + dist.size(), fft_upper);
  const Eigen::Index n_fft = last - dist.data();
  if (n_fft == 0)
    return fit;

  // This resolution kept the worst paired ISE perturbation below 0.7%.
  const size_t num_bins = 256;
  fft::KdeFFT kde_fft(
    dist.head(n_fft), h, 0.0, fft_upper, Eigen::VectorXd(), num_bins);
  const double scale =
    static_cast<double>(n_fft) / static_cast<double>(dist.size());
  const Eigen::VectorXd f0 = scale * kde_fft.kde_drv(0);
  const Eigen::VectorXd f1 = scale * kde_fft.kde_drv(1);
  const Eigen::VectorXd f2 = scale * kde_fft.kde_drv(2);
  for (Eigen::Index j = 0; j < eval_dist.size(); ++j) {
    // c1 and c2 are the first inverse-moment rows for degrees 1 and 2.
    const double a = eval_dist(j) / h;
    const double phi = K0_ * std::exp(-0.5 * a * a);
    const double mu0 = 0.5 * std::erfc(-a * 0.7071067811865475);
    const double mu1 = -phi;
    const double mu2 = mu0 - a * phi;
    const double mu3 = -(a * a + 2.0) * phi;
    const double mu4 = 3.0 * mu0 - (a * a * a + 3.0 * a) * phi;
    const double det1 = mu0 * mu2 - mu1 * mu1;
    const double c10 = mu2 / det1;
    const double c11 = -mu1 / det1;
    const double det2 = mu0 * (mu2 * mu4 - mu3 * mu3) -
                        mu1 * (mu1 * mu4 - mu2 * mu3) +
                        mu2 * (mu1 * mu3 - mu2 * mu2);
    const double c20 = (mu2 * mu4 - mu3 * mu3) / det2;
    const double c21 = (mu2 * mu3 - mu1 * mu4) / det2;
    const double c22 = (mu1 * mu3 - mu2 * mu2) / det2;

    const double position = eval_dist(j) * num_bins / fft_upper;
    const size_t bin = std::min(static_cast<size_t>(position), num_bins - 1);
    const double fraction = position - static_cast<double>(bin);
    auto interpolate = [&](const Eigen::VectorXd& values) {
      return (1.0 - fraction) * values(bin) + fraction * values(bin + 1);
    };
    fit.density(j) =
      0.5 * ((c10 + c20 + c22) * interpolate(f0) -
             h * (c11 + c21) * interpolate(f1) + h * h * c22 * interpolate(f2));
    fit.influence_num(j) =
      K0_ * (c10 + c20) / (2.0 * static_cast<double>(dist.size()) * h);
  }

  // Remove unstable local-polynomial values.
  for (Eigen::Index j = 0; j < fit.density.size(); ++j) {
    if (!(fit.density(j) > 0.0) || !std::isfinite(fit.density(j)) ||
        !std::isfinite(fit.influence_num(j))) {
      fit.density(j) = 0.0;
      fit.influence_num(j) = 0.0;
    }
  }
  return fit;
}

//! Fuses the transformed bulk density @f$\widehat f_B@f$ with endpoint
//! densities @f$\widehat f_L@f$ and @f$\widehat f_U@f$. Each endpoint density
//! averages local-linear and local-quadratic kernels using one degree-2
//! bandwidth, selected from all distances on bounded support and the closest
//! 75% on one-sided support. Smooth weights @f$w_L@f$ and @f$w_U@f$ are
//! functions of the bulk CDF and shrink at rate @f$\min(1/4,n^{-1/2})@f$:
//! @f[
//! \begin{align*}
//! q &= \min(1/4,n^{-1/2}), & z(p) &= \min(1,p/q), \\[2pt]
//! g(p) &= 1-3z(p)^2+2z(p)^3, \\[2pt]
//! w_L &= g(F_B), & w_U &= g(1-F_B), \\[2pt]
//! \widehat f
//!   &= w_L\widehat f_L+(1-w_L-w_U)\widehat f_B+w_U\widehat f_U, \\[2pt]
//! \nu &= w_L\nu_L+(1-w_L-w_U)\nu_B+w_U\nu_U, \\[2pt]
//! \mathrm{influence} &= \nu/\widehat f.
//! \end{align*}
//! @f]
//! Here @f$\nu@f$ is the influence numerator, i.e., the diagonal kernel
//! contribution on the density scale. Only endpoints classified as finite are
//! repaired.
//! @param x observations on the original scale.
//! @param grid increasing evaluation grid on the original scale.
//! @param influences on input, bulk-fit influences on `grid`; on output,
//!   influences of the fused estimate.
inline void
Kde1d::repair_boundaries(const Eigen::VectorXd& x,
                         const Eigen::VectorXd& grid,
                         Eigen::VectorXd& influences)
{
  Eigen::VectorXd dist_lower;
  Eigen::VectorXd dist_upper;
  bool repair_lower = !std::isnan(xmin_);
  bool repair_upper = !std::isnan(xmax_);
  if (repair_lower) {
    dist_lower = x.array() - xmin_;
    repair_lower = is_finite_endpoint(dist_lower);
  }
  if (repair_upper) {
    dist_upper = xmax_ - x.array();
    repair_upper = is_finite_endpoint(dist_upper);
  }
  if (!repair_lower && !repair_upper)
    return;

  // q is the shrinking probability width of each endpoint weight.
  const double q = std::min(0.25, 1.0 / std::sqrt(x.size()));
  auto endpoint_weight = [&](double probability) {
    const double z = std::min(1.0, probability / q);
    return 1.0 - (3.0 * z * z - 2.0 * z * z * z);
  };

  // Only evaluate endpoint fits where their fusion weights are nonzero.
  const Eigen::VectorXd f_bulk = grid_.get_values();
  const Eigen::VectorXd bulk_cdf = grid_.integrate(grid, true);
  Eigen::VectorXd w_l = Eigen::VectorXd::Zero(grid.size());
  Eigen::VectorXd w_u = Eigen::VectorXd::Zero(grid.size());
  for (Eigen::Index j = 0; j < grid.size(); ++j) {
    if (repair_lower)
      w_l(j) = endpoint_weight(bulk_cdf(j));
    if (repair_upper)
      w_u(j) = endpoint_weight(1.0 - bulk_cdf(j));
  }
  const Eigen::Index n_lower = (w_l.array() > 0.0).count();
  const Eigen::Index n_upper = (w_u.array() > 0.0).count();

  // Select one shared degree-2 bandwidth, using the local 75% one-sided rule.
  const double bw_fraction =
    (!std::isnan(xmin_) && !std::isnan(xmax_)) ? 1.0 : 0.75;
  Eigen::VectorXd bw_dist = repair_lower ? dist_lower : dist_upper;
  std::sort(bw_dist.data(), bw_dist.data() + bw_dist.size());
  const Eigen::Index n_bw = std::min<Eigen::Index>(
    bw_dist.size(),
    std::max<Eigen::Index>(
      4, static_cast<Eigen::Index>(std::ceil(bw_fraction * bw_dist.size()))));
  bandwidth::PluginBandwidthSelector selector(bw_dist.head(n_bw));
  const double h = selector.select_bandwidth(2) * multiplier_;

  auto fit_endpoint = [&](const Eigen::VectorXd& dist,
                          const Eigen::VectorXd& eval_dist) {
    Eigen::VectorXd dist_sorted = dist;
    std::sort(dist_sorted.data(), dist_sorted.data() + dist_sorted.size());
    return fit_boundary_component(dist_sorted, eval_dist, h);
  };

  BoundaryComponent lower;
  BoundaryComponent upper;
  repair_lower = repair_lower && n_lower > 0;
  repair_upper = repair_upper && n_upper > 0;
  if (repair_lower) {
    lower = fit_endpoint(dist_lower, grid.head(n_lower).array() - xmin_);
    repair_lower = lower.density.size() > 0;
  }
  if (repair_upper) {
    Eigen::VectorXd upper_eval_dist =
      (xmax_ - grid.tail(n_upper).array()).reverse();
    upper = fit_endpoint(dist_upper, upper_eval_dist);
    repair_upper = upper.density.size() > 0;
  }
  if (!repair_lower && !repair_upper)
    return;

  // Fuse densities and influence numerators with the same weights.
  Eigen::VectorXd f(grid.size());
  Eigen::VectorXd infl_num(grid.size());
  for (Eigen::Index j = 0; j < grid.size(); ++j) {
    const double w_b = 1.0 - w_l(j) - w_u(j);
    f(j) = w_b * f_bulk(j);
    infl_num(j) = w_b * f_bulk(j) * influences(j);
    if (w_l(j) > 0.0) {
      f(j) += w_l(j) * lower.density(j);
      infl_num(j) += w_l(j) * lower.influence_num(j);
    }
    if (w_u(j) > 0.0) {
      const Eigen::Index upper_j = grid.size() - 1 - j;
      f(j) += w_u(j) * upper.density(upper_j);
      infl_num(j) += w_u(j) * upper.influence_num(upper_j);
    }
    influences(j) = f(j) > 0.0 ? infl_num(j) / f(j) : 0.0;
  }
  grid_ = interp::InterpolationGrid(grid, f.cwiseMax(0.0), 3);
}

//  Bandwidth for Kernel Density Estimation
//' @param x vector of observations
//' @param bandwidth bandwidth parameter, NA for automatic selection.
//' @param multiplier bandwidth multiplieriplier.
//' @param discrete whether a jittered estimate is computed.
//' @param weights vector of weights for each observation (can be empty).
//' @param degree polynomial degree.
//' @return the selected bandwidth
//' @noRd
inline double
Kde1d::select_bandwidth(const Eigen::VectorXd& x,
                        double bandwidth,
                        double multiplier,
                        size_t degree,
                        const Eigen::VectorXd& weights) const
{
  if (std::isnan(bandwidth)) {
    bandwidth::PluginBandwidthSelector selector(x, weights);
    bandwidth = selector.select_bandwidth(degree);
  }

  bandwidth *= multiplier;
  if (type_ == VarType::discrete) {
    bandwidth = std::max(bandwidth, 0.5 / 5);
  }

  return bandwidth;
}

inline void
Kde1d::check_xmin_xmax(const double& xmin, const double& xmax) const
{
  if (!std::isnan(xmax) && !std::isnan(xmax) && (xmin > xmax))
    throw std::invalid_argument("xmin must be smaller than xmax");
}

inline void
Kde1d::check_fitted() const
{
  if (grid_.get_grid_points().size() == 0) {
    throw std::runtime_error("You must first fit the KDE to data.");
  }
}

inline void
Kde1d::check_notfitted() const
{
  if (grid_.get_grid_points().size() > 0) {
    throw std::runtime_error(
      "This method can't be used for already fitted objects.");
  }
}

inline void
Kde1d::check_inputs(const Eigen::VectorXd& x,
                    const Eigen::VectorXd& weights) const
{
  if (x.size() == 0)
    throw std::invalid_argument("x must not be empty");

  if ((weights.size() > 0) && (weights.size() != x.size()))
    throw std::invalid_argument("x and weights must have the same size");
}

inline void
Kde1d::check_boundaries(const Eigen::VectorXd& x) const
{
  if ((x.array() < xmin_).any() || (x.array() > xmax_).any()) {
    throw std::invalid_argument("x must be contained in [xmin, xmax].");
  }
}

void
Kde1d::set_interpolation_grid(const interp::InterpolationGrid& grid)
{
  grid_ = grid;
}

void
Kde1d::set_xmin_xmax(double xmin, double xmax)
{
  this->check_notfitted();
  this->check_xmin_xmax(xmin, xmax);
  xmin_ = xmin;
  xmax_ = xmax;
}

std::string
Kde1d::as_str(VarType type) const
{
  std::string type_str;
  switch (type) {
    case VarType::continuous:
      return "continuous";
    case VarType::discrete:
      return "discrete";
    case VarType::zero_inflated:
      return "zero-inflated";
    default:
      throw std::invalid_argument("unknown variable type.");
  }
}

VarType
Kde1d::as_enum(std::string type) const
{
  if ((type == "c") || (type == "cont") || (type == "continuous")) {
    return VarType::continuous;
  } else if ((type == "d") || (type == "disc") || (type == "discrete")) {
    return VarType::discrete;
  } else if ((type == "zi") || (type == "zinfl") || (type == "zero-inflated") ||
             (type == "zero_inflated")) {
    return VarType::zero_inflated;
  } else {
    std::stringstream ss;
    ss << "variable type '" << type << "' unknown; must be one of"
       << "{c, cont, continuous, d, disc, discrete, zi, zinfl, zero-inflated}."
       << std::endl;
    throw std::invalid_argument(ss.str());
  }
  return VarType::continuous;
}

} // end kde1d
