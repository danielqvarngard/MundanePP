#include "magic_reduce/pauli.hpp"

#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void expect(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void test_pauli_1() {
    magic_reduce::Pauli<1> pauli;
    expect(magic_reduce::to_string(pauli) == "I", "Pauli<1> starts as I");
    expect(pauli.weight() == 0, "identity has zero weight");

    pauli.set_x(0);
    expect(pauli.has_x(0) && !pauli.has_z(0), "set_x produces X");
    expect(magic_reduce::to_string(pauli) == "X", "X string rendering");

    pauli.set_y(0);
    expect(pauli.has_x(0) && pauli.has_z(0), "set_y produces Y");

    pauli.set_z(0);
    expect(!pauli.has_x(0) && pauli.has_z(0), "set_z produces Z");

    pauli.clear_qubit(0);
    expect(!pauli.has_x(0) && !pauli.has_z(0), "clear_qubit produces I");
}

void test_pauli_2() {
    magic_reduce::Pauli<2> xx;
    xx.set_x(0);
    xx.set_x(1);

    magic_reduce::Pauli<2> zz;
    zz.set_z(0);
    zz.set_z(1);

    magic_reduce::Pauli<2> xi;
    xi.set_x(0);

    expect(xx.commute(zz), "XX and ZZ commute");
    expect(!xi.commute(zz), "XI and ZZ anticommute");
    expect(xx.weight() == 2, "XX has weight two");
    expect(magic_reduce::to_string(xx) == "XX", "Pauli<2> string rendering");

    auto copy = xx;
    expect(copy.same_pauli(xx), "copied Paulis compare equal");
    copy.flip_z(1);
    expect(!copy.same_pauli(xx), "different Paulis compare unequal");
    copy.flip_z(1);
    expect(copy.same_pauli(xx), "flip_z is reversible");
}

void test_pauli_5() {
    magic_reduce::Pauli<5> pauli;
    pauli.set_x(0);
    pauli.set_y(2);
    pauli.set_z(4);

    expect(pauli.weight() == 3, "Pauli<5> weight spans all labels");
    expect(magic_reduce::to_string(pauli) == "XIYIZ",
           "Pauli<5> string uses ascending qubit order");

    pauli.flip_x(2);
    expect(magic_reduce::to_string(pauli) == "XIZIZ",
           "flip_x changes Y to Z");
    pauli.flip_z(2);
    expect(magic_reduce::to_string(pauli) == "XIIIZ",
           "flip_z changes Z to I");
}

void test_pauli_100() {
    magic_reduce::Pauli<100> pauli;
    pauli.set_x(0);
    pauli.set_y(63);
    pauli.set_z(64);
    pauli.set_x(99);

    expect(magic_reduce::Pauli<100>::blocks == 2,
           "Pauli<100> uses two mask blocks");
    expect(pauli.weight() == 4, "weight crosses a uint64_t boundary");
    expect(pauli.has_x(63) && pauli.has_z(63), "qubit 63 is Y");
    expect(!pauli.has_x(64) && pauli.has_z(64), "qubit 64 is Z");

    magic_reduce::Pauli<100> other;
    other.set_z(0);
    other.set_x(64);
    expect(pauli.commute(other),
           "two symplectic overlaps across blocks commute");

    other.clear_qubit(64);
    expect(!pauli.commute(other), "one symplectic overlap anticommutes");

    auto rendered = magic_reduce::to_string(pauli);
    expect(rendered.size() == 100, "Pauli<100> string has one character per qubit");
    expect(rendered[0] == 'X' && rendered[63] == 'Y' &&
               rendered[64] == 'Z' && rendered[99] == 'X',
           "Pauli<100> string preserves boundary labels");
}

void test_bounds_checking() {
    magic_reduce::Pauli<5> pauli;
    bool threw = false;
    try {
        static_cast<void>(pauli.has_x(5));
    } catch (const std::out_of_range&) {
        threw = true;
    }
    expect(threw, "out-of-range qubit access throws");
}

}  // namespace

int main() {
    try {
        test_pauli_1();
        test_pauli_2();
        test_pauli_5();
        test_pauli_100();
        test_bounds_checking();
    } catch (const std::exception& error) {
        std::cerr << "pauli_tests: " << error.what() << '\n';
        return 1;
    }

    std::cout << "pauli_tests: all tests passed\n";
    return 0;
}
