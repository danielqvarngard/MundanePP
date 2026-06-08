add_library(magic_reduce_headers INTERFACE)
target_compile_features(magic_reduce_headers INTERFACE cxx_std_23)
target_include_directories(magic_reduce_headers INTERFACE
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

add_executable(magic_reduce
    apps/magic_reduce.cpp
)

target_link_libraries(magic_reduce PRIVATE magic_reduce_headers)
