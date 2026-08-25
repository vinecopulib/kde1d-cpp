include_directories(SYSTEM ${external_includes})
add_executable(kde1d-unit-test test.cpp)
target_link_libraries(kde1d-unit-test kde1d)

add_executable(kde1d-invariant-test test_numerical_invariants.cpp)
target_link_libraries(kde1d-invariant-test kde1d)

add_test(NAME kde1d-unit-tests COMMAND $<TARGET_FILE:kde1d-unit-test>)
add_test(NAME kde1d-numerical-invariants COMMAND $<TARGET_FILE:kde1d-invariant-test>)
