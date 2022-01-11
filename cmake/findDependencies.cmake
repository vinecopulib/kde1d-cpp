find_package(Eigen3                       REQUIRED)
find_package(Boost 1.56                   REQUIRED)

set(external_includes ${EIGEN3_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})
