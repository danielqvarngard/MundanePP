# AGENTS.md

## Project: Clifford + Magic Monte Carlo Minimizer

This repository implements a stochastic preconditioner for quantum circuits whose gates are mostly Clifford operations plus arbitrary single-qubit rotations. The first target is a general-purpose **magic minimizer**, not an HQAP-specific solver. HQAP-specific identity-basin heuristics may be added later as optional plugins.

The core goal is:

> Given a circuit, rewrite it into an equivalent or explicitly marked low-loss nearby presentation where the Clifford part absorbs as much structure as possible, while the remaining non-Clifford content is minimized, fused, or localized.

Do not optimize primarily for QASM prettiness or total gate count in the first pass. Optimize for magic structure.

## Preferred architecture

Use a two-language architecture:

1. **Python frontend**

   - Use Qiskit only for QASM import/export and validation.
   - Convert Qiskit circuits into a simple project-owned IR.
   - Convert optimized IR back to QASM.

2. **C++23 core**

   - Implement the Monte Carlo magic minimizer.
   - Use compact bit-packed Pauli representations.
   - Keep the hot loop independent of Qiskit.
   - Support MPI parallel tempering later, but keep a serial executable working first.

## Build assumptions

Use:

- C++23
- CMake
- Python 3
- pytest for Python tests
- Catch2 or simple CTest-based C++ tests
- Optional MPI, behind a build flag

Avoid heavy dependencies in the C++ core unless there is a clear reason.

## Core representation

Represent Pauli strings by X/Z masks. For each qubit q:

| x_q | z_q | Pauli |
| --- | --- | ----- |
| 0   | 0   | I     |
| 1   | 0   | X     |
| 0   | 1   | Z     |
| 1   | 1   | Y     |

For fixed-size kernels, use:

```cpp
template <std::size_t N>
struct Pauli {
    static constexpr std::size_t blocks = (N + 63) / 64;
    std::array<std::uint64_t, blocks> x{};
    std::array<std::uint64_t, blocks> z{};
};
```

A non-Clifford rotation should be represented as a Pauli gadget:

```cpp
template <std::size_t N>
struct MagicGadget {
    Pauli<N> p;
    double theta = 0.0;
    int location = -1;
    int source_gate = -1;
};
```

A rotation around a Pauli string is:

[
R_P(\theta)=\exp(-i\theta P/2).
]

The C++ core should manipulate these gadgets directly.

## Clifford/magic split

The first program should use a pure split:

```text
Clifford skeleton + magic ledger
```

Do not treat T gates as free. A T gate is simply a particular magic rotation, e.g. `Rz(pi/4)`. The first minimizer should reduce and compact magic without forcing the circuit into Clifford+T form.

A later preconditioner may optionally snap to Clifford+T if the downstream solver benefits from that representation.

## First-pass objective

The first Monte Carlo should minimize a magic-oriented energy, not raw gate count.

Useful energy terms:

```text
N_magic:
    number of non-Clifford gadgets

angle_magic_weight:
    distance from Clifford angles, e.g. theta mod pi/2

pauli_support_weight:
    penalty for high-weight Pauli gadgets

magic_spread:
    penalty for magic spread across many spacetime bins

hard_region_count:
    number of arbitrary-angle regions that resist fusion or Cliffordization

fusion_reward:
    reward for bringing compatible gadgets together
```

A simple angle weight is:

```text
d_C(theta) = distance from theta to nearest multiple of pi/2
w(theta)   = sin^2(d_C(theta))
```

Angles equal to multiples of `pi/2` are Clifford and should have zero magic cost. Angles near but not equal to Clifford angles are small residual magic. Angles like `pi/4` are T-type magic and still count as magic.

## Allowed core moves

Implement the Monte Carlo in stages.

### Stage 1 moves

- Move a magic gadget through a neighboring Clifford gate.
- Conjugate the Pauli mask under H, S, Sdg, X, Y, Z, CX, CZ, and SWAP.
- Fuse two gadgets with identical Pauli masks:
  [
  R_P(a)R_P(b)=R_P(a+b).
  ]
- Reduce angles modulo `2*pi`.
- Absorb gadgets whose angle is Clifford, i.e. approximately a multiple of `pi/2`.
- Track and report gadgets whose angle is close to Clifford but not within tolerance.

### Stage 2 moves

- Temperature-dependent Pauli-weight growth.
- Magic compaction into fewer spacetime bins.
- Optional identity-basin scoring, supplied externally.
- Optional approximate local patch scoring.

### Stage 3 moves

- MPI parallel tempering.
- Replica exchange between neighboring temperatures.
- Checkpointing and archive of best states.

Do not start with Stage 3. First make a correct serial minimizer.

## Correctness rules

The first minimizer should preserve circuit semantics exactly unless a move is explicitly marked as approximate. Approximate moves must record their loss and source gates.

For the initial implementation, avoid approximate moves. Implement only exact Clifford conjugation, exact gadget fusion, and exact Clifford-angle absorption.

Every nontrivial algebraic rule must have a unit test.

## C++ Pauli operations to implement and test

Implement:

- `has_x`
- `has_z`
- `flip_x`
- `flip_z`
- `same_pauli`
- `commute`
- `weight`
- `apply_h`
- `apply_s`
- `apply_sdg`
- `apply_x`
- `apply_y`
- `apply_z`
- `apply_cx`
- `apply_cz`
- `apply_swap`

When Clifford conjugation maps `P -> -P`, absorb the sign by flipping the gadget angle:

```text
R_{-P}(theta) = R_P(-theta)
```

Store Pauli masks canonically without a separate negative-P flag unless a later design requires one.

## QASM frontend

Python should:

1. Read QASM with Qiskit.
2. Unroll or normalize to a small basis:

   - h
   - s
   - sdg
   - x
   - y
   - z
   - cx
   - cz
   - swap
   - rx
   - ry
   - rz

3. Emit a simple JSON or binary IR for the C++ core.
4. Read optimized IR from the C++ core.
5. Reconstruct Qiskit `QuantumCircuit`.
6. Export QASM.

Do not use Qiskit inside the C++ Monte Carlo loop.

## Initial repository layout

Use this layout unless there is a clear reason to change it:

```text
.
├── AGENTS.md
├── CMakeLists.txt
├── README.md
├── pyproject.toml
├── src/
│   ├── magic_reduce/
│   │   ├── __init__.py
│   │   ├── qasm_to_ir.py
│   │   ├── ir_to_qasm.py
│   │   └── schema.py
│   └── cpp/
│       ├── CMakeLists.txt
│       ├── include/
│       │   └── magic_reduce/
│       │       ├── pauli.hpp
│       │       ├── gadget.hpp
│       │       ├── gate.hpp
│       │       ├── energy.hpp
│       │       └── monte_carlo.hpp
│       ├── lib/
│       │   ├── pauli.cpp
│       │   ├── energy.cpp
│       │   └── monte_carlo.cpp
│       └── apps/
│           └── magic_reduce.cpp
├── tests/
│   ├── python/
│   └── cpp/
└── examples/
    ├── small.qasm
    └── README.md
```

## CLI expectations

The first C++ executable should support:

```bash
magic_reduce input.ir --sweeps 10000 --temperature 0.1 --seed 123 --output best.ir
```

Eventually add:

```bash
magic_reduce input.ir \
  --temps geometric:0.01:10.0:32 \
  --sweeps 100000 \
  --swap-every 100 \
  --checkpoint-every 1000 \
  --output best.ir
```

MPI support should be optional and introduced only after the serial code is tested.

## Testing priorities

First tests:

1. Pauli mask creation and equality.
2. Pauli commutation via symplectic inner product.
3. H conjugation:

   - X -> Z
   - Z -> X
   - Y -> -Y, with theta sign flip

4. S conjugation:

   - X -> Y
   - Y -> -X
   - Z -> Z

5. CNOT conjugation:

   - X_control -> X_control X_target
   - Z_target -> Z_control Z_target
   - Z_control unchanged
   - X_target unchanged

6. SWAP conjugation swaps Pauli labels on the two qubits.
7. Gadget fusion:

   - same Pauli adds angles
   - angle equal to multiple of pi/2 is marked Clifford

8. A tiny QASM roundtrip through Python IR.

Do not implement MPI until these pass.

## Style

Prefer:

- Small headers with focused responsibilities.
- Value types for Pauli strings and gadgets.
- Explicit names over clever template metaprogramming.
- Deterministic tests with fixed seeds.
- Clear error messages for unsupported gates.

Avoid:

- Qiskit objects in the C++ core.
- Global mutable state.
- Large rewrites without tests.
- Approximate moves mixed silently with exact moves.
- Premature Clifford+T snapping in the first minimizer.

## Scientific intent

The minimizer should be useful as a general-purpose preconditioner. HQAP-specific structure, such as approximate identity centers obtained from MPO scans, should be an optional scoring extension rather than baked into the core algorithm.

The first milestone is not to solve the full challenge. The first milestone is:

> Given a circuit with Clifford gates and arbitrary rotations, produce an equivalent representation with fewer or more localized non-Clifford Pauli rotation gadgets.
