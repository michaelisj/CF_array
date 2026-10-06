#pragma once

#include <cassert>

#include <iostream>
#include <vector>

using MachineWordT = std::uint64_t;
using SizeT = MachineWordT;

#ifdef _DEBUG
#define TRACE(...) printf(__VA_ARGS__)
#else
#define TRACE(...)
#endif


constexpr auto BITS_IN_BYTE = 8;
constexpr auto MACHINE_WORD_BITS = sizeof(MachineWordT) * BITS_IN_BYTE;

constexpr static MachineWordT getLeadingBitsMask(SizeT index) {
  return (1llu << index) - 1;
}
constexpr static MachineWordT getRangeBitsMask(SizeT left, SizeT right) {
  auto mask = ~getLeadingBitsMask(left);
  if (right != 0) {
    mask &= getLeadingBitsMask(right);
  }
  return mask;
}
constexpr MachineWordT divideUp(MachineWordT nom, MachineWordT den) {
  return (nom + den - 1) / den;
}

template <std::uint64_t AlignTo, typename T> constexpr auto align(T value) {
  static_assert((AlignTo & (AlignTo - 1)) == 0, "Must be a natural power of 2");
  return (value + AlignTo - 1) & ~(AlignTo - 1);
}

constexpr SizeT getMSb(MachineWordT word) {
  return MACHINE_WORD_BITS - __builtin_clzll(word);
}

constexpr SizeT lg2(MachineWordT word) { return getMSb(word) - 1; }

constexpr static auto MACHINE_WORD_SHIFT = lg2(MACHINE_WORD_BITS);
constexpr static auto MACHINE_WORD_MASK = (1 << MACHINE_WORD_SHIFT) - 1;
static_assert((1 << MACHINE_WORD_SHIFT) == MACHINE_WORD_BITS);

class Bitvector {

  constexpr static std::pair<SizeT, SizeT> getBitOffsets(SizeT index) {
    return {index >> MACHINE_WORD_SHIFT, index & MACHINE_WORD_MASK};
  }

public:
  Bitvector(std::uint32_t bitCount)
      : m_bits(align<MACHINE_WORD_BITS>(bitCount) / BITS_IN_BYTE)
      , m_bitCount(bitCount) {}

  MachineWordT getWord(SizeT index, SizeT wordSize) const {
    if (0 == wordSize) {
      return 0;
    }

    TRACE("index=%llu wordSize=%llu bitCount=%u\n", index, wordSize,
          m_bitCount);
    assert(index + wordSize <= m_bitCount && wordSize != 0);
    const auto [indexHi, indexLo] = getBitOffsets(index);
    const auto [indexEndHi, indexEndLo] = getBitOffsets(index + wordSize);
    if (indexHi != indexEndHi && indexEndLo != 0) {
      const auto lo =
          (m_bits[indexHi] & ~getLeadingBitsMask(indexLo)) >> indexLo;
      const auto hi = (m_bits[indexEndHi] & getLeadingBitsMask(indexEndLo))
                      << (MACHINE_WORD_BITS - indexLo);
      return lo | hi;
    }

    TRACE("mask=0x%016llx value=0x%016llx %llu\n",
          getRangeBitsMask(indexLo, indexEndLo), m_bits[indexHi], indexEndLo);
    return (m_bits[indexHi] & getRangeBitsMask(indexLo, indexEndLo)) >> indexLo;
  }

  void setWord(SizeT index, SizeT wordSize, MachineWordT value) {
    TRACE("index=%llu wordSize=%llu bitCount=%u\n value=%llu, ", index,
          wordSize, m_bitCount, value);
    assert(index + wordSize <= m_bitCount);
    assert(value < (1 << wordSize));

    if (0 == wordSize) {
      return;
    }

    const auto [indexHi, indexLo] = getBitOffsets(index);
    const auto [indexEndHi, indexEndLo] = getBitOffsets(index + wordSize);
    if (indexHi != indexEndHi && indexEndLo != 0) {
      auto lo = m_bits[indexHi] & getLeadingBitsMask(indexLo);
      lo |= (value << indexLo); // Will erase some bits on the way.
      m_bits[indexHi] = lo;

      auto hi = m_bits[indexEndHi] & ~getLeadingBitsMask(indexEndLo);
      hi |= value >> (MACHINE_WORD_BITS - indexLo);
      m_bits[indexEndHi] = hi;
    } else {
      m_bits[indexHi] =
          (m_bits[indexHi] & ~getRangeBitsMask(indexLo, indexEndLo)) |
          (value << indexLo);
    }
  }

  void setAllBits(bool bit) {
    const std::uint8_t value = -1 * static_cast<std::uint8_t>(bit);
    memset(m_bits.data(), value, m_bits.size() * sizeof(MachineWordT));
  }

private:
  std::vector<MachineWordT> m_bits;
  std::uint32_t m_bitCount;
};
