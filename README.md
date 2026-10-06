# CF Array

A C++ implementation of a cache-friendly CF array structure using a compact bitvector representation with layered storage for values and pointers.

## Overview

The array stores integer values using a two-level structure:
- **Values (B)**: stored across logarithmic layers based on bit-width
- **Pointers (P)**: index into the value layers using iterated-logarithm layering

This achieves space efficiency by encoding values proportional to their magnitude, supporting both linear construction and random updates.

### Random-input CF array

`CF_arr_random.cpp` implements the variant for random inputs:
- **Hierarchy (B_0, ..., B_W)**: a value of height h is stored in the first free cell `B_d[h][i div 2^(h+d)]`, d = 0, ..., W
- **Tree (T)**: values exceeding n, and values rejected at every depth, are stored in a `std::map` (`using Tree`)
- **Parameters (H, D)**: each index keeps the height and depth of its cell, in plain packed arrays

The structure is rebuilt every n/8 replacements to keep the tree small. Both implementations share `Bitvector.h`.

## Build

```sh
g++ -O2 -o a.out CF_arr.cpp -std=c++17
g++ -O2 -o a.out CF_arr_random.cpp -std=c++17
```

For debug output:

```sh
g++ -D_DEBUG -o a.out CF_arr.cpp -std=c++17
```

## Usage

The `main` function demonstrates construction from a vector, reads, and updates with assertions verifying correctness.
