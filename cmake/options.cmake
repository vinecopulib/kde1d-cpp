set(CMAKE_EXPORT_COMPILE_COMMANDS "ON")
set(CMAKE_MACOSX_RPATH 1)

if(CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE Release)
endif()

option(WARNINGS_AS_ERRORS        "Compiler warnings as errors"       "OFF")
option(OPT_ASAN                  "Use adress sanitizer (debug)"      "ON")
option(BUILD_TESTING             "Build tests."                      "ON")
option(BUILD_BENCHMARKS          "Build benchmarks."                 "OFF")
option(CODE_COVERAGE             "Code coverage."                    "OFF")

# `OPT_ASAN` only reaches `CMAKE_CXX_FLAGS_DEBUG`, so it is silently inert in a
# release build, and it cannot ask for anything but the address sanitizer. This
# takes an explicit `-fsanitize=` list and applies it to every build type.
set(SANITIZERS "" CACHE STRING
    "Comma-separated -fsanitize= list, e.g. address,undefined. Empty disables.")
