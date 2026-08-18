// Includes the umbrella header only, the way a downstream project does. If the
// install flattens the headers, this fails to compile; if the exported target
// loses its Eigen dependency, it fails to configure.
#include <kde1d.hpp>

#include <cstdlib>
#include <iostream>

int
main()
{
  Eigen::VectorXd x(7);
  x << 0.1, 0.4, 0.5, 0.9, 1.3, 1.4, 2.0;

  kde1d::Kde1d fit;
  fit.fit(x);

  const Eigen::VectorXd density = fit.pdf(x);
  if (!(density.array() >= 0.0).all()) {
    std::cerr << "consumer: density must be non-negative\n";
    return EXIT_FAILURE;
  }
  std::cout << "consumer ok, bandwidth=" << fit.get_bandwidth() << "\n";
  return EXIT_SUCCESS;
}
