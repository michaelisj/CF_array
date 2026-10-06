#include "Bitvector.h"

class CFArrayValues {
  /**
   * The array of values B with log(n) layers.
   */
public:
  explicit CFArrayValues(std::size_t size)
      : m_bitvector(size * 4), m_layers(lg2(size - 1) + 1) {
    SizeT cursor = 0;
    SizeT cellSize = 0;
    SizeT cells = size;
    while (cells > 0) {
      m_layerStartingPosition.push_back(cursor);
      cursor += cells * cellSize;
      cellSize++;
      cells >>= 1;
    }
    assert(cursor < size * 4);
  }

  MachineWordT getWord(SizeT layer, SizeT index) const {
    auto bitIndex = m_layerStartingPosition[layer] + layer * (index >> layer);
    return m_bitvector.getWord(bitIndex, layer) + (1 << layer);
  }

  void setWord(SizeT layer, SizeT index, MachineWordT value) {
    auto bitIndex = m_layerStartingPosition[layer] + layer * (index >> layer);
    m_bitvector.setWord(bitIndex, layer, value - (1 << layer));
  }

private:
  Bitvector m_bitvector;
  SizeT m_layers;
  std::vector<SizeT> m_layerStartingPosition;
};

class CFArrayPointers {
  /**
   * The array of pointers P with log*(n) layers.
   */
public:
  explicit CFArrayPointers(std::size_t size)
      : m_bitvector(size * 4) // This is an overestimate.
      , m_size(size) {
    SizeT cursor = 0;
    SizeT cellSize = 1;
    SizeT cells = 2 * size;
    const auto highestValue = lg2(size - 1) + 1; // Approximately
    while (true) {
      const auto effectiveCells = cells >> cellSize;
      m_layerIndexAndSize.emplace_back(cursor, cellSize);
      cursor += effectiveCells * cellSize;
      if ((1 << cellSize) > highestValue) {
        break;
      }
      cellSize = (1 << cellSize);
    }
    assert(cursor < size * 4);
    m_bitvector.setAllBits(false);
  }

  MachineWordT get(SizeT index) const {
    if (m_bitvector.getWord(index, 1) == 1) {
      return 0;
    }
    SizeT layer = 1;
    SizeT cellSize = 2;
    while (true) {
      const auto bitIndex = m_layerIndexAndSize[layer].first +
                            cellSize * (index >> (cellSize - 1));
      const auto value = m_bitvector.getWord(bitIndex, cellSize);
      TRACE("index=%llu, bitIndex=%llu values=%llu cellSize=%llu, sz=%zu\n", index,
            bitIndex, value, cellSize, m_layerIndexAndSize.size());
      if (value != 0) {
        return value;
      }
      layer++;
      cellSize = 1 << cellSize;
    }
  }

  void set(SizeT index, MachineWordT value) {
    if (0 == value) {
      m_bitvector.setWord(index, 1, 1);
      return;
    }
    m_bitvector.setWord(index, 1, 0);

    SizeT layer = 1;
    SizeT cellSize = 2;
    while (true) {
      const auto bitIndex = m_layerIndexAndSize[layer].first +
                            cellSize * (index >> (cellSize - 1));
      if (0 == (value >> cellSize)) {
        m_bitvector.setWord(bitIndex, cellSize, value);
        return;
      }
      m_bitvector.setWord(bitIndex, cellSize, 0);
      layer++;
      cellSize = 1 << cellSize;
    }
  }

  void setFast(SizeT layer, SizeT index, MachineWordT value) {
    if (0 == value) {
      m_bitvector.setWord(index, 1, 1);
      return;
    }
    const auto cellSize = m_layerIndexAndSize[layer].second;
    const auto bitIndex =
        m_layerIndexAndSize[layer].first + cellSize * (index >> (cellSize - 1));
    m_bitvector.setWord(bitIndex, cellSize, value);
  }

private:
  Bitvector m_bitvector;
  SizeT m_size;
  std::vector<std::pair<SizeT, SizeT>> m_layerIndexAndSize;
};

class CFArray {
public:
  explicit CFArray(const std::vector<MachineWordT> &values)
      : m_B(values.size()) // This is an overestimate.
      , m_P(values.size()) {
    linearConstruction(values);
  }

  void setAllValues(const std::vector<MachineWordT> &values) {
    for (int i = 0; i < values.size(); ++i) {
      set(i, values[i]);
    }
  }

  void linearConstruction(const std::vector<MachineWordT> &values) {
    const auto layers = lg2(values.size() - 1) + 1;
    SizeT pIndex = 1;
    SizeT pSize = 2;

    m_layerMap.push_back(0);
    for (int i = 1; i < layers; ++i) {
      if (i >= (1 << pSize)) {
        pIndex++;
        pSize = 1 << pSize;
      }
      m_layerMap.push_back(pIndex);
    }

    for (int i = 0; i < values.size(); ++i) {
      setOptimized(i, values[i]);
    }
  }

  MachineWordT get(SizeT index) const {
    const auto layer = m_P.get(index);
    return m_B.getWord(layer, index);
  }

  void set(SizeT index, MachineWordT value) {
    const auto layer = lg2(value);
    m_B.setWord(layer, index, value);
    m_P.set(index, layer);
  }

  void setOptimized(SizeT index, MachineWordT value) {
    const auto bLayer = lg2(value);
    m_B.setWord(bLayer, index, value);
    const auto pLayer = m_layerMap[bLayer];
    m_P.setFast(pLayer, index, bLayer);
  }

private:
  CFArrayValues m_B;
  CFArrayPointers m_P;
  std::vector<MachineWordT> m_layerMap;
};

int main() {
  std::vector<MachineWordT> sample = {
      2, 1, 1, 7,  1, 3, 1, 12, 2, 1, 1, 7,  1, 3, 1, 24, 2, 1, 1, 7,
      1, 3, 1, 12, 2, 1, 1, 7,  1, 3, 1, 48, 2, 1, 1, 7,  1, 3, 1, 12};
  CFArray cfarr(sample);

  for (int i = 0; i < sample.size(); ++i) {
    assert(cfarr.get(i) == sample[i]);
  }
  printf("Read test passed.\n");

  std::vector<std::pair<MachineWordT, MachineWordT>> updates = {{0, 3}, {7, 14}};
  for (const auto &[index, value] : updates) {
    cfarr.set(index, value);
    sample[index] = value;
  }
  for (int i = 0; i < sample.size(); ++i) {
    assert(cfarr.get(i) == sample[i]);
  }
  printf("Update test passed.\n");

  return 0;
}
