# magic-reduce

`magic-reduce` is a stochastic preconditioner for quantum circuits composed
primarily of Clifford gates and arbitrary rotations. Its goal is to produce an
equivalent representation with fewer or more localized non-Clifford Pauli
rotation gadgets.

The project uses a two-language architecture:

- A Python frontend will use Qiskit for QASM import, export, and validation,
  converting circuits to and from a project-owned intermediate representation.
- A C++23 core will implement the performance-sensitive magic minimizer without
  depending on Qiskit.

This repository currently contains only the initial buildable skeleton.

## Build

```bash
cmake -S . -B build
cmake --build build
./build/src/cpp/magic_reduce --help
```

## Python package

The empty Python package can be installed in editable mode with:

```bash
python -m pip install -e .
```
