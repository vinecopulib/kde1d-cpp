add_library(kde1d INTERFACE)
target_include_directories(kde1d INTERFACE
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
        )
# The public headers include both <Eigen/Dense> and
# <boost/math/distributions.hpp>, so both belong on the interface -- a consumer
# of the installed package has no other way to learn about them.
target_link_libraries(kde1d INTERFACE Eigen3::Eigen Boost::headers)

if(BUILD_TESTING)
    set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
    add_subdirectory(test)
endif(BUILD_TESTING)

if(BUILD_BENCHMARKS)
    add_subdirectory(benchmark)
endif(BUILD_BENCHMARKS)

if(BUILD_DIAGNOSTICS)
    add_subdirectory(diagnostic)
endif(BUILD_DIAGNOSTICS)

# Related to exports for linux/mac and code coverage
####
# Installation

# Layout. This works for all platforms:
#   * <prefix>/lib/cmake/kde1d
#   * <prefix>/include/
set(config_install_dir "lib/cmake/${PROJECT_NAME}")
set(include_install_dir "include")

set(generated_dir "${CMAKE_CURRENT_BINARY_DIR}/generated")

# Configuration
set(version_config "${generated_dir}/${PROJECT_NAME}ConfigVersion.cmake")
set(project_config "${generated_dir}/${PROJECT_NAME}Config.cmake")
set(targets_export_name "${PROJECT_NAME}Targets")


# Include module with fuction 'write_basic_package_version_file'
include(CMakePackageConfigHelpers)

# Configure '<PROJECT-NAME>ConfigVersion.cmake'
# Note: PROJECT_VERSION is used as a VERSION
write_basic_package_version_file(
        "${version_config}" COMPATIBILITY SameMajorVersion
)

# Configure '<PROJECT-NAME>Config.cmake'
# Use variables:
#   * targets_export_name
#   * PROJECT_NAME
configure_package_config_file(
        "cmake/templates/Config.cmake.in"
        "${project_config}"
        INSTALL_DESTINATION "${config_install_dir}"
        PATH_VARS include_install_dir
)

# Targets:
install(TARGETS kde1d EXPORT "${targets_export_name}")


# Install the headers as a directory so the layout is preserved. Globbing them
# was wrong twice over: `GLOB_RECURSE .../include/kde1d.hpp` matches both the
# umbrella `include/kde1d.hpp` *and* `include/kde1d/kde1d.hpp`, and installing
# that list flat put the implementation header on top of the umbrella -- so a
# consumer including <kde1d.hpp> got the implementation, whose quoted includes
# then resolved against the wrong directory. A glob also would not notice a
# header added later without re-running CMake.
install(
        DIRECTORY "${PROJECT_SOURCE_DIR}/include/"
        DESTINATION "${include_install_dir}"
        FILES_MATCHING PATTERN "*.hpp"
)


# Config
#   * <prefix>/lib/cmake/kde1d/kde1dConfig.cmake
#   * <prefix>/lib/cmake/kde1d/kde1dConfigVersion.cmake
install(
        FILES "${project_config}" "${version_config}"
        DESTINATION "${config_install_dir}"
)

# Config
#   * <prefix>/lib/cmake/kde1d/kde1dTargets.cmake
install(
        EXPORT "${targets_export_name}"
        DESTINATION "${config_install_dir}"
)

# Install the export set for code coverage
if(NOT WIN32 AND CMAKE_BUILD_TYPE STREQUAL "Debug" AND BUILD_TESTING AND CODE_COVERAGE)
    include(cmake/codeCoverage.cmake)
    file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/coverage)
    setup_target_for_coverage(${PROJECT_NAME}_coverage
                              "ctest --output-on-failure"
                              coverage)
    add_dependencies(${PROJECT_NAME}_coverage
                     kde1d-unit-test
                     kde1d-invariant-test)
endif()
