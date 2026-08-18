# A minimal downstream project, built against the *installed* package. Guards
# the install layout and the exported target's dependencies: neither is
# exercised by building this repo's own tests, which use the source tree.
cmake_minimum_required(VERSION 3.10)

project(kde1d_consumer CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(kde1d REQUIRED)

add_executable(consumer main.cpp)
target_link_libraries(consumer kde1d)
