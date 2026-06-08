#pragma once

#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>

namespace magic_reduce {

template <std::size_t N> class Pauli {
public:
  static_assert(N > 0, "Pauli must contain at least one qubit");

  static constexpr std::size_t num_qubits = N;
  static constexpr std::size_t blocks = (N + 63) / 64;

  [[nodiscard]] bool has_x(std::size_t qubit) const {
    check_qubit(qubit);
    return (x_[block_index(qubit)] & bit_mask(qubit)) != 0;
  }

  [[nodiscard]] bool has_z(std::size_t qubit) const {
    check_qubit(qubit);
    return (z_[block_index(qubit)] & bit_mask(qubit)) != 0;
  }

  void flip_x(std::size_t qubit) {
    check_qubit(qubit);
    x_[block_index(qubit)] ^= bit_mask(qubit);
  }

  void flip_z(std::size_t qubit) {
    check_qubit(qubit);
    z_[block_index(qubit)] ^= bit_mask(qubit);
  }

  // Set the local Pauli on this qubit to exactly X, Z, or Y.
  void set_x(std::size_t qubit) {
    check_qubit(qubit);
    const auto block = block_index(qubit);
    const auto mask = bit_mask(qubit);
    x_[block] |= mask;
    z_[block] &= ~mask;
  }

  void set_z(std::size_t qubit) {
    check_qubit(qubit);
    const auto block = block_index(qubit);
    const auto mask = bit_mask(qubit);
    x_[block] &= ~mask;
    z_[block] |= mask;
  }

  void set_y(std::size_t qubit) {
    check_qubit(qubit);
    const auto block = block_index(qubit);
    const auto mask = bit_mask(qubit);
    x_[block] |= mask;
    z_[block] |= mask;
  }

  void clear_qubit(std::size_t qubit) {
    check_qubit(qubit);
    const auto block = block_index(qubit);
    const auto mask = bit_mask(qubit);
    x_[block] &= ~mask;
    z_[block] &= ~mask;
  }

  void swap_qubits(std::size_t a, std::size_t b) {
    check_qubit(a);
    check_qubit(b);

    if (a == b) {
      return;
    }

    const bool ax = has_x(a);
    const bool az = has_z(a);
    const bool bx = has_x(b);
    const bool bz = has_z(b);

    clear_qubit(a);
    clear_qubit(b);

    set_local_bits(a, bx, bz);
    set_local_bits(b, ax, az);
  }

  [[nodiscard]] std::uint64_t x_block(std::size_t block) const {
    check_block(block);
    return x_[block] & valid_bits_mask(block);
  }

  [[nodiscard]] std::uint64_t z_block(std::size_t block) const {
    check_block(block);
    return z_[block] & valid_bits_mask(block);
  }

  void clear_unused_bits() noexcept {
    const auto mask = valid_bits_mask(blocks - 1);
    x_[blocks - 1] &= mask;
    z_[blocks - 1] &= mask;
  }

  [[nodiscard]] bool operator==(const Pauli &other) const noexcept {
    for (std::size_t block = 0; block < blocks; ++block) {
      if (x_block_unchecked(block) != other.x_block_unchecked(block)) {
        return false;
      }
      if (z_block_unchecked(block) != other.z_block_unchecked(block)) {
        return false;
      }
    }
    return true;
  }

  [[nodiscard]] bool operator!=(const Pauli &other) const noexcept {
    return !(*this == other);
  }

  [[nodiscard]] std::size_t weight() const noexcept {
    std::size_t result = 0;

    for (std::size_t block = 0; block < blocks; ++block) {
      result += static_cast<std::size_t>(
          std::popcount(x_block_unchecked(block) | z_block_unchecked(block)));
    }

    return result;
  }

  [[nodiscard]] bool commute(const Pauli &other) const noexcept {
    unsigned int parity = 0;

    for (std::size_t block = 0; block < blocks; ++block) {
      parity ^= static_cast<unsigned int>(
          std::popcount(x_block_unchecked(block) &
                        other.z_block_unchecked(block)) &
          1U);

      parity ^= static_cast<unsigned int>(
          std::popcount(z_block_unchecked(block) &
                        other.x_block_unchecked(block)) &
          1U);
    }

    return parity == 0;
  }

private:
  std::array<std::uint64_t, blocks> x_{};
  std::array<std::uint64_t, blocks> z_{};

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

  [[nodiscard]] std::uint64_t
  x_block_unchecked(std::size_t block) const noexcept {
    return x_[block] & valid_bits_mask(block);
  }

  [[nodiscard]] std::uint64_t
  z_block_unchecked(std::size_t block) const noexcept {
    return z_[block] & valid_bits_mask(block);
  }

  void set_local_bits(std::size_t qubit, bool x_bit, bool z_bit) {
    const auto block = block_index(qubit);
    const auto mask = bit_mask(qubit);

    if (x_bit) {
      x_[block] |= mask;
    } else {
      x_[block] &= ~mask;
    }

    if (z_bit) {
      z_[block] |= mask;
    } else {
      z_[block] &= ~mask;
    }
  }

  static void check_qubit(std::size_t qubit) {
    if (qubit >= N) {
      throw std::out_of_range("Pauli qubit index out of range");
    }
  }

  static void check_block(std::size_t block) {
    if (block >= blocks) {
      throw std::out_of_range("Pauli block index out of range");
    }
  }
};

template <std::size_t N> struct PauliHash {
  [[nodiscard]] std::size_t operator()(const Pauli<N> &pauli) const noexcept {
    std::size_t seed = 0xcbf29ce484222325ULL;

    for (std::size_t block = 0; block < Pauli<N>::blocks; ++block) {
      hash_combine(seed, pauli.x_block(block));
      hash_combine(seed, pauli.z_block(block));
    }

    return seed;
  }

private:
  static void hash_combine(std::size_t &seed, std::uint64_t value) noexcept {
    // SplitMix64-style mixing. Good enough for hash-table grouping of Pauli
    // masks.
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    value = value ^ (value >> 31);

    seed ^= static_cast<std::size_t>(value) + 0x9e3779b97f4a7c15ULL +
            (seed << 6) + (seed >> 2);
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
