#include "magic_reduce/gadget.hpp"

#include <cmath>
#include <exception>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <string>

namespace {

constexpr double pi = std::numbers::pi_v<double>;
constexpr double theta = 0.375;
constexpr double epsilon = 1e-12;

void expect(bool condition, const std::string &message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

void expect_near(double actual, double expected, const std::string &message) {
  if (std::abs(actual - expected) > epsilon) {
    throw std::runtime_error(message);
  }
}

template <std::size_t N>
void expect_gadget(const magic_reduce::MagicGadget<N> &gadget,
                   const std::string &pauli, double expected_theta,
                   const std::string &message) {
  expect(magic_reduce::to_string(gadget.p) == pauli,
         message + ": Pauli mask");
  expect_near(gadget.theta, expected_theta, message + ": angle");
}

void test_angle_utilities() {
  using magic_reduce::approximately_clifford_angle;
  using magic_reduce::distance_to_clifford_angle;
  using magic_reduce::normalize_angle;

  expect_near(normalize_angle(0.0), 0.0, "zero normalizes to zero");
  expect_near(normalize_angle(pi), pi, "pi remains pi");
  expect_near(normalize_angle(-pi), pi, "-pi normalizes to pi");
  expect_near(normalize_angle(5.0 * pi), pi,
              "angle larger than two pi normalizes");
  expect_near(normalize_angle(-9.0 * pi / 2.0), -pi / 2.0,
              "large negative angle normalizes");

  expect_near(distance_to_clifford_angle(0.0), 0.0,
              "zero is a Clifford angle");
  expect_near(distance_to_clifford_angle(3.0 * pi / 2.0), 0.0,
              "multiple of pi over two has zero distance");
  expect_near(distance_to_clifford_angle(pi / 4.0), pi / 4.0,
              "pi over four is maximally distant");

  const double close = pi / 2.0 + 1e-9;
  expect(approximately_clifford_angle(close, 2e-9),
         "nearby multiple is approximately Clifford");
  expect(!approximately_clifford_angle(close, 5e-10),
         "tight tolerance rejects nearby multiple");
}

void test_h_conjugation() {
  magic_reduce::MagicGadget<1> x;
  x.p.set_x(0);
  x.theta = theta;
  magic_reduce::apply_h(x, 0);
  expect_gadget(x, "Z", theta, "H X H = Z");

  magic_reduce::MagicGadget<1> z;
  z.p.set_z(0);
  z.theta = theta;
  magic_reduce::apply_h(z, 0);
  expect_gadget(z, "X", theta, "H Z H = X");

  magic_reduce::MagicGadget<1> y;
  y.p.set_y(0);
  y.theta = theta;
  magic_reduce::apply_h(y, 0);
  expect_gadget(y, "Y", -theta, "H Y H = -Y");
}

void test_s_conjugation() {
  magic_reduce::MagicGadget<1> x;
  x.p.set_x(0);
  x.theta = theta;
  magic_reduce::apply_s(x, 0);
  expect_gadget(x, "Y", theta, "S X Sdg = Y");

  magic_reduce::MagicGadget<1> y;
  y.p.set_y(0);
  y.theta = theta;
  magic_reduce::apply_s(y, 0);
  expect_gadget(y, "X", -theta, "S Y Sdg = -X");

  magic_reduce::MagicGadget<1> z;
  z.p.set_z(0);
  z.theta = theta;
  magic_reduce::apply_s(z, 0);
  expect_gadget(z, "Z", theta, "S Z Sdg = Z");

  magic_reduce::MagicGadget<1> sdg_x;
  sdg_x.p.set_x(0);
  sdg_x.theta = theta;
  magic_reduce::apply_sdg(sdg_x, 0);
  expect_gadget(sdg_x, "Y", -theta, "Sdg X S = -Y");

  magic_reduce::MagicGadget<1> sdg_y;
  sdg_y.p.set_y(0);
  sdg_y.theta = theta;
  magic_reduce::apply_sdg(sdg_y, 0);
  expect_gadget(sdg_y, "X", theta, "Sdg Y S = X");

  magic_reduce::MagicGadget<1> sdg_z;
  sdg_z.p.set_z(0);
  sdg_z.theta = theta;
  magic_reduce::apply_sdg(sdg_z, 0);
  expect_gadget(sdg_z, "Z", theta, "Sdg Z S = Z");
}

void test_pauli_conjugation() {
  magic_reduce::MagicGadget<1> x_on_z;
  x_on_z.p.set_z(0);
  x_on_z.theta = theta;
  magic_reduce::apply_x(x_on_z, 0);
  expect_gadget(x_on_z, "Z", -theta, "X Z X = -Z");

  magic_reduce::MagicGadget<1> x_on_x;
  x_on_x.p.set_x(0);
  x_on_x.theta = theta;
  magic_reduce::apply_x(x_on_x, 0);
  expect_gadget(x_on_x, "X", theta, "X X X = X");

  magic_reduce::MagicGadget<1> x_on_y;
  x_on_y.p.set_y(0);
  x_on_y.theta = theta;
  magic_reduce::apply_x(x_on_y, 0);
  expect_gadget(x_on_y, "Y", -theta, "X Y X = -Y");

  magic_reduce::MagicGadget<1> y_on_x;
  y_on_x.p.set_x(0);
  y_on_x.theta = theta;
  magic_reduce::apply_y(y_on_x, 0);
  expect_gadget(y_on_x, "X", -theta, "Y X Y = -X");

  magic_reduce::MagicGadget<1> y_on_y;
  y_on_y.p.set_y(0);
  y_on_y.theta = theta;
  magic_reduce::apply_y(y_on_y, 0);
  expect_gadget(y_on_y, "Y", theta, "Y Y Y = Y");

  magic_reduce::MagicGadget<1> y_on_z;
  y_on_z.p.set_z(0);
  y_on_z.theta = theta;
  magic_reduce::apply_y(y_on_z, 0);
  expect_gadget(y_on_z, "Z", -theta, "Y Z Y = -Z");

  magic_reduce::MagicGadget<1> z_on_x;
  z_on_x.p.set_x(0);
  z_on_x.theta = theta;
  magic_reduce::apply_z(z_on_x, 0);
  expect_gadget(z_on_x, "X", -theta, "Z X Z = -X");

  magic_reduce::MagicGadget<1> z_on_y;
  z_on_y.p.set_y(0);
  z_on_y.theta = theta;
  magic_reduce::apply_z(z_on_y, 0);
  expect_gadget(z_on_y, "Y", -theta, "Z Y Z = -Y");

  magic_reduce::MagicGadget<1> z_on_z;
  z_on_z.p.set_z(0);
  z_on_z.theta = theta;
  magic_reduce::apply_z(z_on_z, 0);
  expect_gadget(z_on_z, "Z", theta, "Z Z Z = Z");
}

void test_cx_conjugation() {
  magic_reduce::MagicGadget<2> x_control;
  x_control.p.set_x(0);
  x_control.theta = theta;
  magic_reduce::apply_cx(x_control, 0, 1);
  expect_gadget(x_control, "XX", theta, "CX maps X control to XX");

  magic_reduce::MagicGadget<2> z_target;
  z_target.p.set_z(1);
  z_target.theta = theta;
  magic_reduce::apply_cx(z_target, 0, 1);
  expect_gadget(z_target, "ZZ", theta, "CX maps Z target to ZZ");

  magic_reduce::MagicGadget<2> z_control;
  z_control.p.set_z(0);
  z_control.theta = theta;
  magic_reduce::apply_cx(z_control, 0, 1);
  expect_gadget(z_control, "ZI", theta, "CX leaves Z control unchanged");

  magic_reduce::MagicGadget<2> x_target;
  x_target.p.set_x(1);
  x_target.theta = theta;
  magic_reduce::apply_cx(x_target, 0, 1);
  expect_gadget(x_target, "IX", theta, "CX leaves X target unchanged");

  magic_reduce::MagicGadget<2> mixed;
  mixed.p.set_x(0);
  mixed.p.set_z(1);
  mixed.theta = theta;
  magic_reduce::apply_cx(mixed, 0, 1);
  expect_gadget(mixed, "YY", -theta,
                "CX mixed Pauli triggers canonical sign flip");
}

void test_cz_conjugation() {
  magic_reduce::MagicGadget<2> x_a;
  x_a.p.set_x(0);
  x_a.theta = theta;
  magic_reduce::apply_cz(x_a, 0, 1);
  expect_gadget(x_a, "XZ", theta, "CZ maps X_a to X_a Z_b");

  magic_reduce::MagicGadget<2> x_b;
  x_b.p.set_x(1);
  x_b.theta = theta;
  magic_reduce::apply_cz(x_b, 0, 1);
  expect_gadget(x_b, "ZX", theta, "CZ maps X_b to Z_a X_b");

  magic_reduce::MagicGadget<2> zz;
  zz.p.set_z(0);
  zz.p.set_z(1);
  zz.theta = theta;
  magic_reduce::apply_cz(zz, 0, 1);
  expect_gadget(zz, "ZZ", theta, "CZ leaves Z terms unchanged");

  magic_reduce::MagicGadget<2> mixed;
  mixed.p.set_x(0);
  mixed.p.set_y(1);
  mixed.theta = theta;
  magic_reduce::apply_cz(mixed, 0, 1);
  expect_gadget(mixed, "YX", -theta, "CZ maps X_a Y_b to -Y_a X_b");
}

void test_swap_conjugation() {
  magic_reduce::MagicGadget<2> gadget;
  gadget.p.set_y(0);
  gadget.p.set_z(1);
  gadget.theta = theta;
  gadget.location = 7;
  gadget.source_gate = 11;

  magic_reduce::apply_swap(gadget, 0, 1);

  expect_gadget(gadget, "ZY", theta,
                "SWAP exchanges labels without changing sign");
  expect(gadget.location == 7 && gadget.source_gate == 11,
         "conjugation preserves gadget metadata");
}

} // namespace

int main() {
  try {
    test_angle_utilities();
    test_h_conjugation();
    test_s_conjugation();
    test_pauli_conjugation();
    test_cx_conjugation();
    test_cz_conjugation();
    test_swap_conjugation();
  } catch (const std::exception &error) {
    std::cerr << "gadget_tests: " << error.what() << '\n';
    return 1;
  }

  std::cout << "gadget_tests: all tests passed\n";
  return 0;
}
