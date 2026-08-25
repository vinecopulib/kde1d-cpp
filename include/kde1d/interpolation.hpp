#pragma once

#include "tools.hpp"
#include <Eigen/Dense>
#include <vector>

namespace kde1d {

namespace interp {

//! A class for cubic spline interpolation in one dimension
//!
//! The class is used for implementing kernel estimators. It makes storing the
//! observations obsolete and allows for fast numerical integration.
class InterpolationGrid
{
public:
  InterpolationGrid() {}

  InterpolationGrid(const Eigen::VectorXd& grid_points,
                    const Eigen::VectorXd& values,
                    int norm_times);

  void normalize(int times);

  Eigen::VectorXd interpolate(const Eigen::VectorXd& x) const;

  Eigen::VectorXd integrate(const Eigen::VectorXd& u,
                            bool normalize = false) const;

  Eigen::VectorXd get_values() const { return values_; }
  Eigen::VectorXd get_grid_points() const { return grid_points_; }
  double get_grid_max() const { return grid_points_[grid_points_.size() - 1]; }
  double get_grid_min() const { return grid_points_[0]; }

private:
  // Utility functions for spline Interpolation
  double cubic_poly(const double& x, const Eigen::Vector4d& a) const;
  double cubic_indef_integral(const double& x, const Eigen::Vector4d& a) const;
  double cubic_integral(const double& lower,
                        const double& upper,
                        const Eigen::Vector4d& a) const;
  size_t binary_search(const double& x) const;
  size_t find_cell(const double& x) const;
  Eigen::Vector4d find_cell_coefs(const size_t& k) const;
  void update_cell_coefs();
  void update_cumulative_integrals();
  void update_cell_lookup();

  Eigen::VectorXd grid_points_;
  Eigen::VectorXd values_;
  Eigen::Matrix<double, 4, Eigen::Dynamic> cell_coefs_;
  Eigen::VectorXd cumulative_integrals_;
  std::vector<size_t> cell_lookup_;
};

//! Constructor
//!
//! @param grid_points an ascending sequence of grid points.
//! @param values a vector of values of same length as grid_points.
//! @param norm_times how many times the normalization routine should run.
inline InterpolationGrid::InterpolationGrid(const Eigen::VectorXd& grid_points,
                                            const Eigen::VectorXd& values,
                                            int norm_times)
{
  if (grid_points.size() != values.size())
    throw std::invalid_argument(
      "grid_points and values must be of equal length");

  grid_points_ = grid_points;
  values_ = values;
  this->normalize(norm_times);
  update_cell_lookup();
}

//! renormalizes the estimate to integrate to one
//!
//! @param times how many times the normalization routine should run.
inline void
InterpolationGrid::normalize(int times)
{
  for (int k = 0; k < times; ++k) {
    double integral = 0.0;
    for (Eigen::Index cell = 0; cell < grid_points_.size() - 1; ++cell) {
      integral += cubic_integral(0.0, 1.0, find_cell_coefs(cell)) *
                  (grid_points_(cell + 1) - grid_points_(cell));
    }
    values_ /= integral;
  }
  update_cell_coefs();
  update_cumulative_integrals();
}

//! Interpolation
//! @param x vector of evaluation points.
inline Eigen::VectorXd
InterpolationGrid::interpolate(const Eigen::VectorXd& x) const
{
  auto interpolate_one = [&](const double& xx) {
    size_t k = find_cell(xx);
    double xev =
      (xx - grid_points_(k)) / (grid_points_(k + 1) - grid_points_(k));

    // use Gaussian tail for extrapolation
    if (xev <= 0) {
      return values_(k) * std::exp(-0.5 * xev * xev);
    } else if (xev >= 1) {
      return values_(k + 1) * std::exp(-0.5 * (xev - 1) * (xev - 1));
    }

    return cubic_poly(xev, Eigen::Vector4d(cell_coefs_.col(k)));
  };

  return tools::unaryExpr_or_nan(x, interpolate_one);
}

//! Integration along the grid
//!
//! @param x a vector  of evaluation points
//! @param normalize whether to normalize the integral to a maximum value of 1.
inline Eigen::VectorXd
InterpolationGrid::integrate(const Eigen::VectorXd& x, bool normalize) const
{
  Eigen::VectorXd res(x.size());
  for (Eigen::Index i = 0; i < x.size(); ++i) {
    if (std::isnan(x(i))) {
      res(i) = x(i);
    } else if (x(i) <= grid_points_(0)) {
      res(i) = 0.0;
    } else if (x(i) >= grid_points_(grid_points_.size() - 1)) {
      res(i) = cumulative_integrals_(grid_points_.size() - 1);
    } else {
      const size_t cell = find_cell(x(i));
      const double cell_width = grid_points_(cell + 1) - grid_points_(cell);
      const double position = (x(i) - grid_points_(cell)) / cell_width;
      res(i) = cumulative_integrals_(cell) +
               cubic_integral(
                 0.0, position, Eigen::Vector4d(cell_coefs_.col(cell))) *
                 cell_width;
    }
  }

  return normalize ? res / cumulative_integrals_(grid_points_.size() - 1)
                   : res;
}

// ---------------- Utility functions for spline interpolation ----------------

//! Evaluate a cubic polynomial
//!
//! @param x evaluation point.
//! @param a polynomial coefficients
inline double
InterpolationGrid::cubic_poly(const double& x, const Eigen::Vector4d& a) const
{
  double x2 = x * x;
  double x3 = x2 * x;
  return a(0) + a(1) * x + a(2) * x2 + a(3) * x3;
}

//! Indefinite integral of a cubic polynomial
//!
//! @param x evaluation point.
//! @param a polynomial coefficients.
inline double
InterpolationGrid::cubic_indef_integral(const double& x,
                                        const Eigen::Vector4d& a) const
{
  double x2 = x * x;
  double x3 = x2 * x;
  double x4 = x3 * x;
  return a(0) * x + a(1) / 2.0 * x2 + a(2) / 3.0 * x3 + a(3) / 4.0 * x4;
}

//! Definite integral of a cubic polynomial
//!
//! @param lower lower limit of the integral.
//! @param upper upper limit of the integral.
//! @param a polynomial coefficients.
inline double
InterpolationGrid::cubic_integral(const double& lower,
                                  const double& upper,
                                  const Eigen::Vector4d& a) const
{
  return cubic_indef_integral(upper, a) - cubic_indef_integral(lower, a);
}

inline size_t
InterpolationGrid::binary_search(const double& x) const
{
  size_t low = 0, high = grid_points_.size() - 1;
  size_t mid;
  while (low < high - 1) {
    mid = low + (high - low) / 2;
    if (x < grid_points_(mid))
      high = mid;
    else
      low = mid;
  }

  return low;
}

inline size_t
InterpolationGrid::find_cell(const double& x) const
{
  const double relative_position =
    (x - grid_points_(0)) /
    (grid_points_(grid_points_.size() - 1) - grid_points_(0));
  const double bounded_position =
    std::min(std::max(relative_position, 0.0), 1.0);
  const size_t bucket = std::min(
    static_cast<size_t>(bounded_position *
                        static_cast<double>(cell_lookup_.size())),
    cell_lookup_.size() - 1);
  size_t cell = cell_lookup_[bucket];
  while ((cell < static_cast<size_t>(grid_points_.size() - 2)) &&
         (grid_points_(cell + 1) <= x)) {
    ++cell;
  }
  return cell;
}

//! Calculate coefficients for cubic intrpolation spline
//!
//! @param k the cell index.
inline Eigen::Vector4d
InterpolationGrid::find_cell_coefs(const size_t& k) const
{
  // indices for cell and neighboring grid points
  long int k0 =
    std::max(static_cast<long int>(k) - 1, static_cast<long int>(0));
  long k2 = k + 1;
  long k3 = std::min(static_cast<long int>(k + 2),
                     static_cast<long int>(grid_points_.size() - 1));

  double dt0 = grid_points_(k) - grid_points_(k0);
  double dt1 = grid_points_(k2) - grid_points_(k);
  double dt2 = grid_points_(k3) - grid_points_(k2);

  // compute tangents when parameterized in (t1,t2)
  // for smooth extrapolation, derivative is set to zero at boundary
  double dx1 = 0.0, dx2 = 0.0;
  if (dt0 > 0) {
    dx1 = (values_(k) - values_(k0)) / dt0;
    dx1 -= (values_(k2) - values_(k0)) / (dt0 + dt1);
    dx1 += (values_(k2) - values_(k)) / dt1;
  }
  if (dt2 > 0) {
    dx2 = (values_(k2) - values_(k)) / dt1;
    dx2 -= (values_(k3) - values_(k)) / (dt1 + dt2);
    dx2 += (values_(k3) - values_(k2)) / dt2;
  }

  // rescale tangents for parametrization in (0,1)
  dx1 *= dt1;
  dx2 *= dt1;

  // ensure positivity (Schmidt and Hess, DOI:10.1007/bf01934097)
  dx1 = std::max(dx1, -3 * values_(k));
  dx2 = std::min(dx2, 3 * values_(k2));

  // compute coefficents
  Eigen::Vector4d a;
  a(0) = values_(k);
  a(1) = dx1;
  a(2) = -3 * (values_(k) - values_(k2)) - 2 * dx1 - dx2;
  a(3) = 2 * (values_(k) - values_(k2)) + dx1 + dx2;

  return a;
}

inline void
InterpolationGrid::update_cell_coefs()
{
  cell_coefs_.resize(4, grid_points_.size() - 1);
  for (Eigen::Index k = 0; k < grid_points_.size() - 1; ++k)
    cell_coefs_.col(k) = find_cell_coefs(k);
}

inline void
InterpolationGrid::update_cumulative_integrals()
{
  cumulative_integrals_ = Eigen::VectorXd::Zero(grid_points_.size());
  for (Eigen::Index cell = 0; cell < grid_points_.size() - 1; ++cell) {
    cumulative_integrals_(cell + 1) =
      cumulative_integrals_(cell) +
      cubic_integral(
        0.0, 1.0, Eigen::Vector4d(cell_coefs_.col(cell))) *
        (grid_points_(cell + 1) - grid_points_(cell));
  }
}

inline void
InterpolationGrid::update_cell_lookup()
{
  cell_lookup_.resize(128);
  const double grid_range =
    grid_points_(grid_points_.size() - 1) - grid_points_(0);
  for (size_t bucket = 0; bucket < cell_lookup_.size(); ++bucket) {
    cell_lookup_[bucket] = binary_search(
      grid_points_(0) + grid_range * static_cast<double>(bucket) /
                          static_cast<double>(cell_lookup_.size()));
  }
}

} // end kde1d::interp

} // end kde1d
