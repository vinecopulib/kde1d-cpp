find_package(Eigen3                       REQUIRED)
find_package(Boost 1.56                   REQUIRED)

set(external_includes ${EIGEN3_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})

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
