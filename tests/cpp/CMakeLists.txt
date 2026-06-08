add_executable(pauli_tests
    test_pauli.cpp
)

target_link_libraries(pauli_tests PRIVATE magic_reduce_headers)

add_test(NAME pauli_tests COMMAND pauli_tests)
