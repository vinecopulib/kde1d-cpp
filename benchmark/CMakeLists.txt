add_executable(kde1d-benchmark benchmark.cpp)
target_include_directories(kde1d-benchmark PRIVATE ${PROJECT_SOURCE_DIR}/test)
target_link_libraries(kde1d-benchmark kde1d)
set_target_properties(
  kde1d-benchmark
  PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
