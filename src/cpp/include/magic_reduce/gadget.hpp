#pragma once

#include "magic_reduce/pauli.hpp"

#include <cmath>
#include <cstddef>
#include <numbers>
#include <stdexcept>
#include <string>

namespace magic_reduce {

// A gadget represents R_P(theta) = exp(-i theta P / 2).
//
// The conjugation functions below implement forward Clifford conjugation:
// P -> G P G^\dagger, so G R_P(theta) G^\dagger =
// R_{G P G^\dagger}(theta). If the result is -P', the Pauli mask remains
// canonical and theta is negated because R_{-P'}(theta) = R_{P'}(-theta).
template <std::size_t N> struct MagicGadget {
  Pauli<N> p;
  double theta = 0.0;
  int location = -1;
  int source_gate = -1;
};

// Returns the equivalent angle in (-pi, pi]. In particular, both pi and -pi
// normalize to pi. Non-finite inputs follow the standard library remainder
// behavior and produce NaN.
[[nodiscard]] inline double normalize_angle(double theta) noexcept {
  constexpr double pi = std::numbers::pi_v<double>;
  constexpr double two_pi = 2.0 * pi;

  double normalized = std::remainder(theta, two_pi);
  if (normalized <= -pi) {
    normalized = pi;
  }
  return normalized;
}

template <std::size_t N>
[[nodiscard]] MagicGadget<N>
make_magic_gadget(Pauli<N> p, double theta, int location = -1,
                  int source_gate = -1) {
  return MagicGadget<N>{p, normalize_angle(theta), location, source_gate};
}

// Returns the shortest angular distance to any integer multiple of pi / 2.
[[nodiscard]] inline double
distance_to_clifford_angle(double theta) noexcept {
  constexpr double half_pi = std::numbers::pi_v<double> / 2.0;
  return std::abs(std::remainder(theta, half_pi));
}

[[nodiscard]] inline bool approximately_clifford_angle(double theta,
                                                        double tolerance) {
  if (tolerance < 0.0) {
    throw std::invalid_argument("angle tolerance must be non-negative");
  }
  return distance_to_clifford_angle(theta) <= tolerance;
}

namespace detail {

template <std::size_t N>
void flip_theta_sign(MagicGadget<N> &gadget) noexcept {
  gadget.theta = normalize_angle(-gadget.theta);
}

template <std::size_t N>
void set_local_pauli(Pauli<N> &pauli, std::size_t qubit, bool x, bool z) {
  if (x && z) {
    pauli.set_y(qubit);
  } else if (x) {
    pauli.set_x(qubit);
  } else if (z) {
    pauli.set_z(qubit);
  } else {
    pauli.clear_qubit(qubit);
  }
}

inline void check_distinct_qubits(std::size_t a, std::size_t b,
                                  const char *gate_name) {
  if (a == b) {
    throw std::invalid_argument(std::string{gate_name} +
                                " requires distinct qubits");
  }
}

} // namespace detail

template <std::size_t N>
void apply_h(MagicGadget<N> &gadget, std::size_t qubit) {
  const bool x = gadget.p.has_x(qubit);
  const bool z = gadget.p.has_z(qubit);

  if (x && z) {
    detail::flip_theta_sign(gadget);
  }
  detail::set_local_pauli(gadget.p, qubit, z, x);
}

template <std::size_t N>
void apply_s(MagicGadget<N> &gadget, std::size_t qubit) {
  const bool x = gadget.p.has_x(qubit);
  const bool z = gadget.p.has_z(qubit);

  if (x && z) {
    detail::flip_theta_sign(gadget);
  }
  detail::set_local_pauli(gadget.p, qubit, x, z ^ x);
}

template <std::size_t N>
void apply_sdg(MagicGadget<N> &gadget, std::size_t qubit) {
  const bool x = gadget.p.has_x(qubit);
  const bool z = gadget.p.has_z(qubit);

  if (x && !z) {
    detail::flip_theta_sign(gadget);
  }
  detail::set_local_pauli(gadget.p, qubit, x, z ^ x);
}

template <std::size_t N>
void apply_x(MagicGadget<N> &gadget, std::size_t qubit) {
  if (gadget.p.has_z(qubit)) {
    detail::flip_theta_sign(gadget);
  }
}

template <std::size_t N>
void apply_y(MagicGadget<N> &gadget, std::size_t qubit) {
  if (gadget.p.has_x(qubit) != gadget.p.has_z(qubit)) {
    detail::flip_theta_sign(gadget);
  }
}

template <std::size_t N>
void apply_z(MagicGadget<N> &gadget, std::size_t qubit) {
  if (gadget.p.has_x(qubit)) {
    detail::flip_theta_sign(gadget);
  }
}

template <std::size_t N>
void apply_cx(MagicGadget<N> &gadget, std::size_t control,
              std::size_t target) {
  detail::check_distinct_qubits(control, target, "CX");

  const bool x_c = gadget.p.has_x(control);
  const bool z_c = gadget.p.has_z(control);
  const bool x_t = gadget.p.has_x(target);
  const bool z_t = gadget.p.has_z(target);

  if (x_c && z_t && (x_t == z_c)) {
    detail::flip_theta_sign(gadget);
  }
  if (x_c) {
    gadget.p.flip_x(target);
  }
  if (z_t) {
    gadget.p.flip_z(control);
  }
}

template <std::size_t N>
void apply_cz(MagicGadget<N> &gadget, std::size_t a, std::size_t b) {
  detail::check_distinct_qubits(a, b, "CZ");

  const bool x_a = gadget.p.has_x(a);
  const bool z_a = gadget.p.has_z(a);
  const bool x_b = gadget.p.has_x(b);
  const bool z_b = gadget.p.has_z(b);

  if (x_a && x_b && (z_a ^ z_b)) {
    detail::flip_theta_sign(gadget);
  }
  if (x_b) {
    gadget.p.flip_z(a);
  }
  if (x_a) {
    gadget.p.flip_z(b);
  }
}

template <std::size_t N>
void apply_swap(MagicGadget<N> &gadget, std::size_t a, std::size_t b) {
  detail::check_distinct_qubits(a, b, "SWAP");
  gadget.p.swap_qubits(a, b);
}

} // namespace magic_reduce
