#pragma once

#include "dpik.hpp"
#include "interpolation.hpp"
#include "stats.hpp"
#include "tools.hpp"
#include <cmath>
#include <functional>

namespace kde1d {

//! Local-polynomial density estimation in 1-d.
class Kde1d
{
public:
  // constructors
  Kde1d(double bw = NAN,
        double mult = 1.0,
        double xmin = NAN,
        double xmax = NAN,
        size_t deg = 2);

  Kde1d(const interp::InterpolationGrid& grid,
        double xmin = NAN,
        double xmax = NAN);

  void fit(const Eigen::VectorXd& x,
           const Eigen::VectorXd& weights = Eigen::VectorXd());

  void fit(const Eigen::VectorXi& x,
           const Eigen::VectorXd& weights = Eigen::VectorXd());

  // statistical functions
  Eigen::VectorXd pdf(const Eigen::VectorXd& x) const;
  Eigen::VectorXd pdf(const Eigen::VectorXi& x) const;
  Eigen::VectorXd cdf(const Eigen::VectorXd& x) const;
  Eigen::VectorXd cdf(const Eigen::VectorXi& x) const;
  Eigen::VectorXd quantile(const Eigen::VectorXd& x) const;
  Eigen::VectorXd simulate(size_t n, const std::vector<int>& seeds = {}) const;

  // getters
  Eigen::VectorXd get_values() const { return grid_.get_values(); }
  Eigen::VectorXd get_grid_points() const { return grid_.get_grid_points(); }
  double get_bw() const { return bw_; }
  double get_mult() const { return mult_; }
  double get_deg() const { return deg_; }
  double get_xmin() const { return xmin_; }
  double get_xmax() const { return xmax_; }
  double get_edf() const { return edf_; }
  double get_loglik() const { return loglik_; }

  std::string str() const
  {
    std::stringstream ss;
    ss << "Kde1d("
       << "bw=" << bw_ << ", mult=" << mult_ << ", xmin=" << xmin_
       << ", xmax=" << xmax_ << ", deg=" << deg_ << ")";
    return ss.str();
  }

private:
  // data members
  interp::InterpolationGrid grid_;
  double xmin_;
  double xmax_;
  bool is_discrete_{ false };
  double bw_{ NAN };
  double mult_{ 1.0 };
  size_t deg_{ 2 };
  double loglik_{ NAN };
  double edf_{ NAN };
  static constexpr double K0_ = 0.3989425;

  // private methods
  void fit_internal(const Eigen::VectorXd& x,
                    const Eigen::VectorXd& weights = Eigen::VectorXd());
  Eigen::VectorXd pdf_continuous(const Eigen::VectorXd& x) const;
  Eigen::VectorXd cdf_continuous(const Eigen::VectorXd& x) const;
  Eigen::VectorXd quantile_continuous(const Eigen::VectorXd& x) const;
  Eigen::VectorXd pdf_discrete(const Eigen::VectorXi& x) const;
  Eigen::VectorXd cdf_discrete(const Eigen::VectorXi& x) const;
  Eigen::VectorXd quantile_discrete(const Eigen::VectorXi& x) const;

  Eigen::VectorXd kern_gauss(const Eigen::VectorXd& x);
  Eigen::MatrixXd fit_lp(const Eigen::VectorXd& x,
                         const Eigen::VectorXd& grid,
                         const Eigen::VectorXd& weights);
  double calculate_infl(const size_t& n,
                        const double& f0,
                        const double& b,
                        const double& bw,
                        const double& s,
                        const double& weight);
  Eigen::VectorXd boundary_transform(const Eigen::VectorXd& x,
                                     bool inverse = false);
  Eigen::VectorXd boundary_correct(const Eigen::VectorXd& x,
                                   const Eigen::VectorXd& fhat);
  Eigen::VectorXd construct_grid_points(const Eigen::VectorXd& x);
  Eigen::VectorXd finalize_grid(Eigen::VectorXd& grid_points);
  double select_bw(const Eigen::VectorXd& x,
                   double bw,
                   double mult,
                   size_t deg,
                   const Eigen::VectorXd& weights) const;
};

//! constructor for fitting the density estimate.
//! @param bw positive bandwidth parameter (NaN means automatic selection).
//! @param mult bandwidth multiplier.
//! @param xmin lower bound for the support of the density, `NaN` means no
//!   boundary.
//! @param xmax upper bound for the support of the density, `NaN` means no
//!   boundary.
//! @param deg order of the local polynomial.
inline Kde1d::Kde1d(double bw,
                    double mult,
                    double xmin,
                    double xmax,
                    size_t deg)
  : xmin_(xmin)
  , xmax_(xmax)
  , bw_(bw)
  , mult_(mult)
  , deg_(deg)
{
  if (deg_ > 2)
    throw std::runtime_error("deg must not be larger than 2.");
}

//! @param x vector of observations
//! @param weights vector of weights for each observation (optional).
inline void
Kde1d::fit(const Eigen::VectorXd& x, const Eigen::VectorXd& weights)
{
  is_discrete_ = false;
  fit_internal(x, weights);
}

//! @param x vector of observations
//! @param weights vector of weights for each observation (optional).
inline void
Kde1d::fit(const Eigen::VectorXi& x, const Eigen::VectorXd& weights)
{
  is_discrete_ = true;
  xmin_ = NAN;
  xmax_ = NAN;
  fit_internal(x.cast<double>(), weights);
}

//! @param x vector of observations
//! @param weights vector of weights for each observation (optional).
inline void
Kde1d::fit_internal(const Eigen::VectorXd& x, const Eigen::VectorXd& weights)
{
  if (weights.size() > 0 && (weights.size() != x.size()))
    throw std::runtime_error("x and weights must have the same size.");

  // preprocessing for nans and jittering
  Eigen::VectorXd xx = x;
  Eigen::VectorXd w = weights;
  tools::remove_nans(xx, w);
  if (w.size() > 0)
    w /= w.mean();
  if (is_discrete_)
    xx = stats::equi_jitter(xx);
  xx = boundary_transform(xx);

  // bandwidth selection
  bw_ = select_bw(xx, bw_, mult_, deg_, w);

  // fit model and evaluate in transformed domain
  Eigen::VectorXd grid_points = construct_grid_points(xx);
  Eigen::MatrixXd fitted = fit_lp(xx, boundary_transform(grid_points), w);

  // correct estimated density for transformation
  Eigen::VectorXd values = boundary_correct(grid_points, fitted.col(0));

  // move boundary points to xmin/xmax
  grid_points = finalize_grid(grid_points);

  // construct interpolation grid
  // (3 iterations for normalization to a proper density)
  grid_ = interp::InterpolationGrid(grid_points, values, 3);

  // calculate log-likelihood of final estimate
  loglik_ = grid_.interpolate(x).cwiseMax(1e-20).array().log().sum();

  // calculate effective degrees of freedom
  interp::InterpolationGrid infl_grid(
    grid_points, fitted.col(1).cwiseMin(2.0).cwiseMax(0), 0);
  edf_ = infl_grid.interpolate(x).sum();
}

//! construct model from an already fit interpolation grid.
//! @param grid the interpolation grid.
//! @param xmin lower bound for the support of the density, `NaN` means no
//!   boundary.
//! @param xmax upper bound for the support of the density, `NaN` means no
//!   boundary.
inline Kde1d::Kde1d(const interp::InterpolationGrid& grid,
                    double xmin,
                    double xmax)
  : grid_(grid)
  , xmin_(xmin)
  , xmax_(xmax)
{}

//! computes the pdf of the kernel density estimate by interpolation.
//! @param x vector of evaluation points.
//! @return a vector of pdf values.
inline Eigen::VectorXd
Kde1d::pdf(const Eigen::VectorXd& x) const
{
  return !is_discrete_ ? pdf_continuous(x) : pdf_discrete(x.cast<int>());
}

//! computes the pdf of the kernel density estimate by interpolation.
//! @param x vector of evaluation points.
//! @return a vector of pdf values.
inline Eigen::VectorXd
Kde1d::pdf(const Eigen::VectorXi& x) const
{
  return !is_discrete_ ? pdf_continuous(x.cast<double>()) : pdf_discrete(x);
}

inline Eigen::VectorXd
Kde1d::pdf_continuous(const Eigen::VectorXd& x) const
{
  Eigen::VectorXd fhat = grid_.interpolate(x);
  if (!std::isnan(xmin_)) {
    fhat = (x.array() < xmin_).select(Eigen::VectorXd::Zero(x.size()), fhat);
  }
  if (!std::isnan(xmax_)) {
    fhat = (x.array() > xmax_).select(Eigen::VectorXd::Zero(x.size()), fhat);
  }

  auto trunc = [](const double& xx) { return std::max(xx, 0.0); };
  return tools::unaryExpr_or_nan(fhat, trunc);
  ;
}

inline Eigen::VectorXd
Kde1d::pdf_discrete(const Eigen::VectorXi& x) const
{
  auto fhat = pdf_continuous(x.cast<double>());
  // normalize
  auto mx = std::lround(grid_.get_grid_max());
  auto mn = std::lround(grid_.get_grid_min());
  Eigen::VectorXd lvs = Eigen::VectorXd::LinSpaced(mx - mn + 1, mn, mx);
  fhat /= grid_.interpolate(lvs).sum();

  return fhat;
}

//! computes the cdf of the kernel density estimate by numerical
//! integration.
//! @param x vector of evaluation points.
//! @return a vector of cdf values.
inline Eigen::VectorXd
Kde1d::cdf(const Eigen::VectorXd& x) const
{
  return !is_discrete_ ? cdf_continuous(x) : cdf_discrete(x.cast<int>());
}

//! computes the cdf of the kernel density estimate by numerical
//! integration.
//! @param x vector of evaluation points.
//! @return a vector of cdf values.
inline Eigen::VectorXd
Kde1d::cdf(const Eigen::VectorXi& x) const
{
  return !is_discrete_ ? cdf_continuous(x.cast<double>()) : cdf_discrete(x);
}

inline Eigen::VectorXd
Kde1d::cdf_continuous(const Eigen::VectorXd& x) const
{
  return grid_.integrate(x, /* normalize */ true);
}

inline Eigen::VectorXd
Kde1d::cdf_discrete(const Eigen::VectorXi& x) const
{
  auto mx = std::lround(grid_.get_grid_max());
  auto mn = std::lround(grid_.get_grid_min());
  Eigen::VectorXi lvs = Eigen::VectorXi::LinSpaced(mx - mn + 1, mn, mx);
  auto f_cum = pdf_discrete(lvs);
  for (long i = 1; i <= mx; ++i)
    f_cum(i) += f_cum(i - 1);

  return tools::unaryExpr_or_nan_int(x, [&f_cum](const double& xx) {
    return std::min(1.0, std::max(f_cum(static_cast<size_t>(xx)), 0.0));
  });
}

//! computes the cdf of the kernel density estimate by numerical inversion.
//! @param x vector of evaluation points.
//! @return a vector of quantiles.
inline Eigen::VectorXd
Kde1d::quantile(const Eigen::VectorXd& x) const
{
  if ((x.minCoeff() < 0) | (x.maxCoeff() > 1))
    throw std::runtime_error("probabilities must lie in (0, 1).");
  return !is_discrete_ ? quantile_continuous(x)
                       : quantile_discrete(x.cast<int>());
}

inline Eigen::VectorXd
Kde1d::quantile_continuous(const Eigen::VectorXd& x) const
{
  auto cdf = [&](const Eigen::VectorXd& xx) { return grid_.integrate(xx); };
  auto q =
    tools::invert_f(x, cdf, grid_.get_grid_min(), grid_.get_grid_max(), 35);

  // replace with NaN where the input was NaN
  for (long i = 0; i < x.size(); i++) {
    if (std::isnan(x(i)))
      q(i) = x(i);
  }

  return q;
}

inline Eigen::VectorXd
Kde1d::quantile_discrete(const Eigen::VectorXi& x) const
{
  auto mx = std::lround(grid_.get_grid_max());
  Eigen::VectorXi lvs = Eigen::VectorXi::LinSpaced(mx + 1, 0, mx);
  auto p = cdf_discrete(lvs);
  auto quan = [&](const double& pp) {
    long lv = 0;
    while ((pp >= p(lv)) && (lv < mx))
      lv++;
    return lvs(lv);
  };

  Eigen::VectorXd out(x.size());
  for (long i = 0; i < x.size(); ++i) {
    if (std::isnan(x(i))) {
      out(i) = std::numeric_limits<double>::quiet_NaN();
    } else {
      out(i) = quan(x(i));
    }
  }

  return out;
}

//! simulates data from the model.
//! @param n the number of observations to simulate.
//! @param seeds an optional vector of seeds.
//! @return simulated observations from the kernel density.
inline Eigen::VectorXd
Kde1d::simulate(size_t n, const std::vector<int>& seeds) const
{
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
  fft::KdeFFT kde_fft(x, bw_, grid_points(0), grid_points(m - 1), weights);
  Eigen::VectorXd f0 = kde_fft.kde_drv(0);

  Eigen::VectorXd wbin = Eigen::VectorXd::Ones(m);
  if (weights.size()) {
    // compute the average weight per cell
    auto wcount = kde_fft.get_bin_counts();
    auto count = tools::linbin(x,
                               grid_points(0),
                               grid_points(m - 1),
                               m - 1,
                               Eigen::VectorXd::Ones(x.size()));
    wbin = wcount.cwiseQuotient(count);
  }

  Eigen::MatrixXd res(f0.size(), 2);
  res.col(0) = f0;
  res.col(1) = K0_ / (x.size() * bw_) * wbin.cwiseQuotient(f0);
  if (deg_ == 0)
    return res;

  // deg > 0
  Eigen::VectorXd f1 = kde_fft.kde_drv(1);
  Eigen::VectorXd S = Eigen::VectorXd::Constant(f0.size(), bw_);
  Eigen::VectorXd b = f1.cwiseQuotient(f0);
  if (deg_ == 2) {
    Eigen::VectorXd f2 = kde_fft.kde_drv(2);
    // D/R is notation from Hjort and Jones' AoS paper
    Eigen::VectorXd D = f2.cwiseQuotient(f0) - b.cwiseProduct(b);
    Eigen::VectorXd R = 1 / (1.0 + bw_ * bw_ * D.array()).sqrt();
    // this is our notation
    S = (R / bw_).array().pow(2);
    b *= bw_ * bw_;
    res.col(0) = bw_ * S.cwiseSqrt().cwiseProduct(res.col(0));
  }
  res.col(0) = res.col(0).array() * (-0.5 * b.array().pow(2) * S.array()).exp();

  for (size_t k = 0; k < m; k++) {
    // TODO: weights
    res(k, 1) = calculate_infl(x.size(), f0(k), b(k), bw_, S(k), wbin(k));
    if (std::isnan(res(k, 0)))
      res.row(k).setZero();
  }

  return res;
}

//! calculate influence for data point for density estimate based on
//! quantities pre-computed in `fit_lp()`.
inline double
Kde1d::calculate_infl(const size_t& n,
                      const double& f0,
                      const double& b,
                      const double& bw,
                      const double& s,
                      const double& weight)
{
  double M_inverse00;
  double bw2 = std::pow(bw, 2);
  double b2 = std::pow(b, 2);
  if (deg_ == 0) {
    M_inverse00 = 1 / f0;
  } else if (deg_ == 1) {
    Eigen::Matrix2d M;
    M(0, 0) = f0;
    M(0, 1) = bw2 * b * f0;
    M(1, 0) = M(0, 1);
    M(1, 1) = f0 * bw2 + f0 * bw2 * bw2 * b2;
    M_inverse00 = M.inverse()(0, 0);
  } else {
    Eigen::Matrix3d M;
    M(0, 0) = f0;
    M(0, 1) = f0 * b;
    M(1, 0) = M(0, 1);
    M(1, 1) = f0 * bw2 + f0 * b2;
    M(1, 2) = 0.5 * f0 * (3.0 / s * b + b * b2);
    M(2, 1) = M(1, 2);
    M(2, 2) = 0.25 * f0;
    M(2, 2) *= 3.0 / std::pow(s, 2) + 6.0 / s * b2 + b2 * b2;
    M(0, 2) = M(2, 2);
    M(2, 0) = M(2, 2);
    M_inverse00 = M.inverse()(0, 0);
  }

  return K0_ * weight / (n * bw) * M_inverse00;
}

//! transformations for density estimates with bounded support.
//! @param x evaluation points.
//! @param inverse whether the inverse transformation should be applied.
//! @return the transformed evaluation points.
inline Eigen::VectorXd
Kde1d::boundary_transform(const Eigen::VectorXd& x, bool inverse)
{
  Eigen::VectorXd x_new = x;
  if (!inverse) {
    if (!std::isnan(xmin_) & !std::isnan(xmax_)) {
      // two boundaries -> probit transform
      auto rng = xmax_ - xmin_;
      x_new = (x.array() - xmin_ + 5e-5 * rng) / (xmax_ - xmin_ + 1e-4 * rng);
      x_new = stats::qnorm(x_new);
    } else if (!std::isnan(xmin_)) {
      // left boundary -> log transform
      x_new = (1e-5 + x.array() - xmin_).log();
    } else if (!std::isnan(xmax_)) {
      // right boundary -> negative log transform
      x_new = (1e-5 + xmax_ - x.array()).log();
    } else {
      // no boundary -> no transform
    }
  } else {
    if (!std::isnan(xmin_) & !std::isnan(xmax_)) {
      // two boundaries -> probit transform
      auto rng = xmax_ - xmin_;
      x_new = stats::pnorm(x).array() + xmin_ - 5e-5 * rng;
      x_new *= (xmax_ - xmin_ + 1e-4 * rng);
    } else if (!std::isnan(xmin_)) {
      // left boundary -> log transform
      x_new = x.array().exp() + xmin_ - 1e-5;
    } else if (!std::isnan(xmax_)) {
      // right boundary -> negative log transform
      x_new = -(x.array().exp() - xmax_ - 1e-5);
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
  Eigen::VectorXd corr_term(fhat.size());
  if (!std::isnan(xmin_) & !std::isnan(xmax_)) {
    // two boundaries -> probit transform
    auto rng = xmax_ - xmin_;
    corr_term = (x.array() - xmin_ + 5e-5 * rng) / (xmax_ - xmin_ + 1e-4 * rng);
    corr_term = stats::dnorm(stats::qnorm(corr_term));
    corr_term /= (xmax_ - xmin_ + 1e-4 * rng);
    corr_term = 1.0 / corr_term.array().max(1e-6);
  } else if (!std::isnan(xmin_)) {
    // left boundary -> log transform
    corr_term = 1.0 / (1e-5 + x.array() - xmin_).max(1e-6);
  } else if (!std::isnan(xmax_)) {
    // right boundary -> negative log transform
    corr_term = 1.0 / (1e-5 + xmax_ - x.array()).max(1e-6);
  } else {
    // no boundary -> no transform
    corr_term.fill(1.0);
  }

  Eigen::VectorXd f_corr = fhat.cwiseProduct(corr_term);
  if (std::isnan(xmin_) & !std::isnan(xmax_))
    f_corr.reverseInPlace();

  return f_corr;
}

//! constructs a grid later used for interpolation
//! @param x vector of observations.
//! @return a grid of size 50.
inline Eigen::VectorXd
Kde1d::construct_grid_points(const Eigen::VectorXd& x)
{
  Eigen::VectorXd rng(2);
  rng << x.minCoeff(), x.maxCoeff();
  if (std::isnan(xmin_) && std::isnan(xmax_)) {
    rng(0) -= 4 * bw_;
    rng(1) += 4 * bw_;
  }
  auto zgrid = Eigen::VectorXd::LinSpaced(401, rng(0), rng(1));
  return boundary_transform(zgrid, true);
}

//! moves the boundary points of the grid to xmin/xmax (if non-NaN).
//! @param grid_points the grid points.
inline Eigen::VectorXd
Kde1d::finalize_grid(Eigen::VectorXd& grid_points)
{
  if (std::isnan(xmin_) & !std::isnan(xmax_))
    grid_points.reverseInPlace();
  if (!std::isnan(xmin_))
    grid_points(0) = xmin_;
  if (!std::isnan(xmax_))
    grid_points(grid_points.size() - 1) = xmax_;

  return grid_points;
}

//  Bandwidth for Kernel Density Estimation
//' @param x vector of observations
//' @param bw bandwidth parameter, NA for automatic selection.
//' @param mult bandwidth multiplier.
//' @param discrete whether a jittered estimate is computed.
//' @param weights vector of weights for each observation (can be empty).
//' @param deg polynomial degree.
//' @return the selected bandwidth
//' @noRd
inline double
Kde1d::select_bw(const Eigen::VectorXd& x,
                 double bw,
                 double mult,
                 size_t deg,
                 const Eigen::VectorXd& weights) const
{
  if (std::isnan(bw)) {
    bw::PluginBandwidthSelector selector(x, weights);
    bw = selector.select_bw(deg);
  }

  bw *= mult;
  if (is_discrete_) {
    bw = std::max(bw, 0.5 / 5);
  }

  return bw;
}

} // end kde1d
