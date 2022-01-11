include_directories(SYSTEM ${external_includes})
add_executable(test test.cpp)
target_link_libraries(test kde1d)
