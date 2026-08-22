add_executable(kde1d-boundary-diagnostic boundary_diagnostic.cpp)
target_link_libraries(kde1d-boundary-diagnostic kde1d)
set_target_properties(
  kde1d-boundary-diagnostic
  PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin
)
