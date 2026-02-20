# CF Array

A C++ implementation of a cache-friendly CF array structure using a compact bitvector representation with layered storage for values and pointers.

## Overview

The array stores integer values using a two-level structure:
- **Values (B)**: stored across logarithmic layers based on bit-width
- **Pointers (P)**: index into the value layers using iterated-logarithm layering

This achieves space efficiency by encoding values proportional to their magnitude, supporting both linear construction and random updates.

## Build

```sh
g++ -O2 -o a.out CF_arr.cpp -std=c++17
```

For debug output:

```sh
g++ -D_DEBUG -o a.out CF_arr.cpp -std=c++17
```

## Usage

The `main` function demonstrates construction from a vector, reads, and updates with assertions verifying correctness.
