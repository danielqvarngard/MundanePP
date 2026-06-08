#include "magic_reduce/pauli.hpp"

#include <cstdint>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_set>

namespace {

void expect(bool condition, const std::string &message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

template <std::size_t N>
void expect_equal_hash(const magic_reduce::Pauli<N> &left,
                       const magic_reduce::Pauli<N> &right,
                       const std::string &message) {
  expect(left == right, message + ": values compare equal");
  expect(magic_reduce::PauliHash<N>{}(left) ==
             magic_reduce::PauliHash<N>{}(right),
         message + ": equal values have equal hashes");
}

void test_pauli_1() {
  magic_reduce::Pauli<1> pauli;
  const magic_reduce::Pauli<1> identity;

  expect(!pauli.has_x(0), "Pauli<1> identity has no X bit");
  expect(!pauli.has_z(0), "Pauli<1> identity has no Z bit");
  expect(pauli.weight() == 0, "Pauli<1> identity has zero weight");
  expect(magic_reduce::to_string(pauli) == "I",
         "Pauli<1> identity renders as I");
  expect_equal_hash(pauli, identity, "Pauli<1> identity");

  pauli.flip_x(0);
  expect(pauli.has_x(0) && !pauli.has_z(0), "flip_x changes I to X");
  expect(pauli.weight() == 1, "X has weight one");
  expect(magic_reduce::to_string(pauli) == "X", "X renders correctly");
  expect(pauli != identity, "X differs from identity");

  pauli.flip_z(0);
  expect(pauli.has_x(0) && pauli.has_z(0), "flip_z changes X to Y");
  expect(magic_reduce::to_string(pauli) == "Y", "Y renders correctly");

  pauli.set_z(0);
  expect(!pauli.has_x(0) && pauli.has_z(0), "set_z replaces Y with Z");

  pauli.set_y(0);
  expect(pauli.has_x(0) && pauli.has_z(0), "set_y replaces Z with Y");

  pauli.set_x(0);
  expect(pauli.has_x(0) && !pauli.has_z(0), "set_x replaces Y with X");

  pauli.clear_qubit(0);
  expect_equal_hash(pauli, identity, "cleared Pauli<1>");

  magic_reduce::Pauli<1> z;
  z.set_z(0);
  magic_reduce::Pauli<1> x;
  x.set_x(0);
  expect(!x.commute(z), "single-qubit X and Z anticommute");
  expect(x.commute(x), "a Pauli commutes with itself");
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
  expect(magic_reduce::to_string(xx) == "XX",
         "Pauli<2> renders in increasing qubit order");

  auto copy = xx;
  expect_equal_hash(copy, xx, "copied Pauli<2>");
  copy.flip_z(1);
  expect(copy != xx, "flipping one mask bit changes equality");
  expect(!(copy == xx), "operator== rejects different masks");
  copy.flip_z(1);
  expect_equal_hash(copy, xx, "reversibly flipped Pauli<2>");

  std::unordered_set<magic_reduce::Pauli<2>, magic_reduce::PauliHash<2>>
      values;
  values.insert(xx);
  values.insert(copy);
  values.insert(zz);
  expect(values.size() == 2,
         "PauliHash and equality deduplicate equal Pauli<2> values");
}

void test_pauli_5() {
  magic_reduce::Pauli<5> pauli;
  pauli.set_x(0);
  pauli.set_y(2);
  pauli.set_z(4);

  expect(pauli.weight() == 3, "Pauli<5> weight counts non-identity labels");
  expect(magic_reduce::to_string(pauli) == "XIYIZ",
         "Pauli<5> string uses increasing qubit order");

  pauli.flip_x(2);
  expect(magic_reduce::to_string(pauli) == "XIZIZ",
         "flip_x changes Y to Z");
  pauli.flip_z(2);
  expect(magic_reduce::to_string(pauli) == "XIIIZ",
         "flip_z changes Z to I");

  magic_reduce::Pauli<5> same;
  same.set_x(0);
  same.set_z(4);
  expect_equal_hash(pauli, same, "independently constructed Pauli<5>");

  same.set_y(3);
  expect(same != pauli, "operator!= detects a Pauli<5> difference");
}

void test_pauli_100() {
  using Pauli100 = magic_reduce::Pauli<100>;

  Pauli100 pauli;
  pauli.set_x(0);
  pauli.set_y(63);
  pauli.set_z(64);
  pauli.set_x(99);

  expect(Pauli100::blocks == 2, "Pauli<100> uses two mask blocks");
  expect(pauli.weight() == 4, "weight crosses a uint64_t boundary");
  expect(pauli.has_x(63) && pauli.has_z(63), "qubit 63 is Y");
  expect(!pauli.has_x(64) && pauli.has_z(64), "qubit 64 is Z");

  const std::uint64_t high_bit = std::uint64_t{1} << 63;
  const std::uint64_t last_valid_bit = std::uint64_t{1} << 35;
  const std::uint64_t final_block_mask =
      (std::uint64_t{1} << 36) - std::uint64_t{1};

  expect(pauli.x_block(0) == (std::uint64_t{1} | high_bit),
         "x_block exposes the complete first X block");
  expect(pauli.z_block(0) == high_bit,
         "z_block exposes the complete first Z block");
  expect(pauli.x_block(1) == last_valid_bit,
         "x_block maps qubit 99 to bit 35 of the final block");
  expect(pauli.z_block(1) == std::uint64_t{1},
         "z_block maps qubit 64 to bit zero of the final block");
  expect((pauli.x_block(1) & ~final_block_mask) == 0 &&
             (pauli.z_block(1) & ~final_block_mask) == 0,
         "public block accessors mask unused Pauli<100> bits");

  Pauli100 other;
  other.set_z(0);
  other.set_x(64);
  expect(pauli.commute(other),
         "two symplectic overlaps across blocks commute");
  other.clear_qubit(64);
  expect(!pauli.commute(other), "one symplectic overlap anticommutes");

  Pauli100 copy = pauli;
  expect_equal_hash(pauli, copy, "copied Pauli<100>");
  copy.clear_qubit(99);
  expect(copy != pauli, "operator!= detects a final-block difference");

  const auto rendered = magic_reduce::to_string(pauli);
  expect(rendered.size() == 100,
         "Pauli<100> string has one character per qubit");
  expect(rendered[0] == 'X' && rendered[63] == 'Y' &&
             rendered[64] == 'Z' && rendered[99] == 'X',
         "Pauli<100> string preserves boundary labels");

  std::unordered_set<Pauli100, magic_reduce::PauliHash<100>> values;
  values.insert(pauli);
  values.insert(Pauli100{pauli});
  values.insert(copy);
  expect(values.size() == 2,
         "PauliHash and equality deduplicate equal Pauli<100> values");
}

void test_bounds_checking() {
  magic_reduce::Pauli<5> pauli;

  bool qubit_threw = false;
  try {
    static_cast<void>(pauli.has_x(5));
  } catch (const std::out_of_range &) {
    qubit_threw = true;
  }
  expect(qubit_threw, "out-of-range qubit access throws");

  bool block_threw = false;
  try {
    static_cast<void>(pauli.x_block(1));
  } catch (const std::out_of_range &) {
    block_threw = true;
  }
  expect(block_threw, "out-of-range block access throws");
}

} // namespace

int main() {
  try {
    test_pauli_1();
    test_pauli_2();
    test_pauli_5();
    test_pauli_100();
    test_bounds_checking();
  } catch (const std::exception &error) {
    std::cerr << "pauli_tests: " << error.what() << '\n';
    return 1;
  }

  std::cout << "pauli_tests: all tests passed\n";
  return 0;
}
