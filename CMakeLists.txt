cmake_minimum_required(VERSION 3.2)

set(CMAKE_CXX_STANDARD 17)

project(kde1d VERSION 1.0)

include(cmake/findDependencies.cmake      REQUIRED)

include(cmake/options.cmake               REQUIRED)

include(cmake/compilerDefOpt.cmake        REQUIRED)

include(cmake/buildTargets.cmake          REQUIRED)

include(cmake/printInfo.cmake             REQUIRED)
