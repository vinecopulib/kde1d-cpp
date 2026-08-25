#include <kde1d.hpp>
#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

namespace {

struct Scenario
{
  std::string name;
  Eigen::VectorXd observations;
  double xmin;
  double xmax;
  double scale;
};

Eigen::VectorXd
probabilities(Eigen::Index size)
{
  return Eigen::VectorXd::LinSpaced(
    size, 0.5 / static_cast<double>(size),
    1.0 - 0.5 / static_cast<double>(size));
}

Scenario
make_scenario(const std::string& name, Eigen::Index size, double scale)
{
  Eigen::VectorXd probability = probabilities(size);
  if (name == "uniform") {
    return { name, probability * scale, 0.0, scale, scale };
  } else if (name == "beta-2-1") {
    return { name, probability.array().sqrt() * scale, 0.0, scale, scale };
  } else if (name == "beta-1-2") {
    return { name,
             (1.0 - (1.0 - probability.array()).sqrt()) * scale,
             0.0,
             scale,
             scale };
  } else if (name == "exponential-left") {
    return { name,
             (-(1.0 - probability.array()).log()) * scale,
             0.0,
             NAN,
             scale };
  }
  return { name,
           (1.0 - probability.array()).log() * scale,
           NAN,
           0.0,
           scale };
}

double
true_pdf(const std::string& name, double position, double scale)
{
  if (name == "uniform")
    return 1.0 / scale;
  if (name == "beta-2-1")
    return 2.0 * position / scale;
  if (name == "beta-1-2")
    return 2.0 * (1.0 - position) / scale;
  return std::exp(-position) / scale;
}

double
true_cdf(const std::string& name, double position)
{
  if (name == "uniform")
    return position;
  if (name == "beta-2-1")
    return position * position;
  if (name == "beta-1-2")
    return 1.0 - (1.0 - position) * (1.0 - position);
  if (name == "exponential-left")
    return 1.0 - std::exp(-position);
  return std::exp(-position);
}

double
evaluation_point(const Scenario& scenario, double position)
{
  if (scenario.name == "exponential-right")
    return -position * scenario.scale;
  return position * scenario.scale;
}

void
run_scenario(const Scenario& scenario,
             size_t degree,
             size_t grid_size,
             double bandwidth,
             const Eigen::VectorXd& positions)
{
  kde1d::Kde1d fit(scenario.xmin,
                   scenario.xmax,
                   "continuous",
                   1.0,
                   bandwidth,
                   degree,
                   grid_size);
  fit.fit(scenario.observations);

  for (Eigen::Index i = 0; i < positions.size(); ++i) {
    const double point = evaluation_point(scenario, positions(i));
    const double estimated_pdf =
      fit.pdf(Eigen::VectorXd::Constant(1, point))(0);
    const double estimated_cdf =
      fit.cdf(Eigen::VectorXd::Constant(1, point))(0);
    std::cout << scenario.name << ',' << scenario.scale << ',' << degree << ','
              << grid_size << ','
              << (std::isnan(bandwidth) ? "automatic"
                                        : std::to_string(bandwidth))
              << ',' << positions(i) << ',' << point << ',' << estimated_pdf
              << ',' << true_pdf(scenario.name, positions(i), scenario.scale)
              << ',' << estimated_cdf << ','
              << true_cdf(scenario.name, positions(i)) << '\n';
  }
}

} // namespace

int
main()
{
  std::cout << std::setprecision(17);
  std::cout << "scenario,scale,degree,grid_size,bandwidth,position,x,"
               "estimated_pdf,true_pdf,estimated_cdf,true_cdf\n";

  Eigen::VectorXd positions(12);
  positions << 0.0, 1e-12, 1e-10, 1e-8, 1e-6, 1e-5, 1e-4, 1e-3,
    1e-2, 1e-1, 0.5, 1.0;
  const std::array<std::string, 5> scenario_names = {
    "uniform", "beta-2-1", "beta-1-2", "exponential-left",
    "exponential-right"
  };

  for (const std::string& scenario_name : scenario_names) {
    for (double scale : { 1e-4, 1.0, 1e4 }) {
      Scenario scenario = make_scenario(scenario_name, 2000, scale);
      for (size_t degree = 0; degree < 3; ++degree) {
        for (size_t grid_size : { 100, 400, 1000 }) {
          for (double bandwidth :
               { 0.15, 0.3, std::numeric_limits<double>::quiet_NaN() }) {
            run_scenario(
              scenario, degree, grid_size, bandwidth, positions);
          }
        }
      }
    }
  }
}
