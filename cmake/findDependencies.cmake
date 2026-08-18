# Eigen in config mode, so it is consumed through its imported target rather
# than an include path. `EIGEN3_INCLUDE_DIR` is set by Eigen's `FindEigen3`
# module but not by its `Eigen3Config.cmake`, so reading that variable silently
# produced an empty include path against an installed Eigen -- point
# `CMAKE_PREFIX_PATH` at the install prefix and `Eigen3::Eigen` carries the
# headers itself.
find_package(Eigen3 3.3                   REQUIRED NO_MODULE)
find_package(Boost 1.56                   REQUIRED)

set(external_includes ${Boost_INCLUDE_DIRS})

# Find doxygen and configure if found
find_package(Doxygen QUIET)
if(DOXYGEN_FOUND)
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in
        ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY
    )
    add_custom_target(doc
        ${DOXYGEN_EXECUTABLE}
        ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
endif (DOXYGEN_FOUND)
