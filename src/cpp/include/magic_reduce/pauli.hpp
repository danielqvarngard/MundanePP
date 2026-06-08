#pragma once

#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>

namespace magic_reduce {

template <std::size_t N> struct Pauli {
  static_assert(N > 0, "Pauli must contain at least one qubit");

  static constexpr std::size_t num_qubits = N;
  static constexpr std::size_t blocks = (N + 63) / 64;

  // Invariant: bits q >= N in the final block should be zero.
  // Most methods mask unused bits defensively, but hash functions should
  // either canonicalize first or mask the final block explicitly.
  std::array<std::uint64_t, blocks> x{};
  std::array<std::uint64_t, blocks> z{};

  [[nodiscard]] bool has_x(std::size_t qubit) const {
    check_qubit(qubit);
    return (x[block_index(qubit)] & bit_mask(qubit)) != 0;
  }

  [[nodiscard]] bool has_z(std::size_t qubit) const {
    check_qubit(qubit);
    return (z[block_index(qubit)] & bit_mask(qubit)) != 0;
  }

  void flip_x(std::size_t qubit) {
    check_qubit(qubit);
    x[block_index(qubit)] ^= bit_mask(qubit);
  }

  void flip_z(std::size_t qubit) {
    check_qubit(qubit);
    z[block_index(qubit)] ^= bit_mask(qubit);
  }

  // Set the local Pauli on this qubit to exactly X, Z, or Y.
  void set_x(std::size_t qubit) {
    check_qubit(qubit);
    const auto block = block_index(qubit);
    const auto mask = bit_mask(qubit);
    x[block] |= mask;
    z[block] &= ~mask;
  }

  void set_z(std::size_t qubit) {
    check_qubit(qubit);
    const auto block = block_index(qubit);
    const auto mask = bit_mask(qubit);
    x[block] &= ~mask;
    z[block] |= mask;
  }

  void set_y(std::size_t qubit) {
    check_qubit(qubit);
    const auto block = block_index(qubit);
    const auto mask = bit_mask(qubit);
    x[block] |= mask;
    z[block] |= mask;
  }

  void clear_qubit(std::size_t qubit) {
    check_qubit(qubit);
    const auto block = block_index(qubit);
    const auto mask = bit_mask(qubit);
    x[block] &= ~mask;
    z[block] &= ~mask;
  }

  void clear_unused_bits() noexcept {
    const auto mask = valid_bits_mask(blocks - 1);
    x[blocks - 1] &= mask;
    z[blocks - 1] &= mask;
  }

  [[nodiscard]] bool same_pauli(const Pauli &other) const noexcept {
    for (std::size_t block = 0; block < blocks; ++block) {
      const auto mask = valid_bits_mask(block);
      if (((x[block] ^ other.x[block]) & mask) != 0 ||
          ((z[block] ^ other.z[block]) & mask) != 0) {
        return false;
      }
    }
    return true;
  }

  [[nodiscard]] bool operator==(const Pauli &other) const noexcept {
    return same_pauli(other);
  }

  [[nodiscard]] bool operator!=(const Pauli &other) const noexcept {
    return !same_pauli(other);
  }

  [[nodiscard]] std::size_t weight() const noexcept {
    std::size_t result = 0;
    for (std::size_t block = 0; block < blocks; ++block) {
      result += static_cast<std::size_t>(
          std::popcount((x[block] | z[block]) & valid_bits_mask(block)));
    }
    return result;
  }

  [[nodiscard]] bool commute(const Pauli &other) const noexcept {
    unsigned int parity = 0;

    for (std::size_t block = 0; block < blocks; ++block) {
      const auto mask = valid_bits_mask(block);

      parity ^= static_cast<unsigned int>(
          std::popcount((x[block] & other.z[block]) & mask) & 1U);

      parity ^= static_cast<unsigned int>(
          std::popcount((z[block] & other.x[block]) & mask) & 1U);
    }

    return parity == 0;
  }

private:
  static constexpr std::size_t block_index(std::size_t qubit) noexcept {
    return qubit / 64;
  }

  static constexpr std::uint64_t bit_mask(std::size_t qubit) noexcept {
    return std::uint64_t{1} << (qubit % 64);
  }

  static constexpr std::uint64_t valid_bits_mask(std::size_t block) noexcept {
    if constexpr (N % 64 == 0) {
      return ~std::uint64_t{0};
    }

    if (block + 1 == blocks) {
      return (std::uint64_t{1} << (N % 64)) - 1;
    }

    return ~std::uint64_t{0};
  }

  static void check_qubit(std::size_t qubit) {
    if (qubit >= N) {
      throw std::out_of_range("Pauli qubit index out of range");
    }
  }
};

// Returns characters in increasing qubit-index order: q0, q1, ..., qN-1.
template <std::size_t N>
[[nodiscard]] std::string to_string(const Pauli<N> &pauli) {
  std::string result;
  result.reserve(N);

  for (std::size_t qubit = 0; qubit < N; ++qubit) {
    const bool has_x = pauli.has_x(qubit);
    const bool has_z = pauli.has_z(qubit);

    if (has_x && has_z) {
      result.push_back('Y');
    } else if (has_x) {
      result.push_back('X');
    } else if (has_z) {
      result.push_back('Z');
    } else {
      result.push_back('I');
    }
  }

  return result;
}

} // namespace magic_reduce
