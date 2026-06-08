cmake_minimum_required(VERSION 3.20)

project(
    magic-reduce
    VERSION 0.1.0
    DESCRIPTION "Clifford and magic Monte Carlo minimizer"
    LANGUAGES CXX
)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(CTest)

add_subdirectory(src/cpp)

if(BUILD_TESTING)
    add_subdirectory(tests/cpp)
endif()
