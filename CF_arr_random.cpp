#include "Bitvector.h"

#include <algorithm>
#include <map>
#include <random>

using Tree = std::map<SizeT, MachineWordT>;

class CFArrayHierarchy {
  /**
   * The arrays B_0, ..., B_W, where W + 1 = layers. In B_d, layer h has
   * ceil(n / 2^(h+d)) cells of h + 1 bits. A cell keeps a height-h value with
   * its leading bit, which serves as the occupancy bit: a free cell reads as 0.
   */
public:
  CFArrayHierarchy(SizeT size, SizeT layers)
      : m_layers(layers)
      , m_rowStartingPosition(getRowStartingPositions(size, layers))
      , m_bitvector(m_rowStartingPosition.back()) {}

  MachineWordT getWord(SizeT depth, SizeT height, SizeT index) const {
    return m_bitvector.getWord(getBitIndex(depth, height, index), height + 1);
  }

  void setWord(SizeT depth, SizeT height, SizeT index, MachineWordT value) {
    m_bitvector.setWord(getBitIndex(depth, height, index), height + 1, value);
  }

private:
  /**
   * The bit offsets of the rows B_d[h], followed by the total size.
   */
  static std::vector<SizeT> getRowStartingPositions(SizeT size, SizeT layers) {
    std::vector<SizeT> positions = {0};
    for (SizeT depth = 0; depth < layers; ++depth) {
      for (SizeT height = 0; height < layers; ++height) {
        const auto cells = divideUp(size, 1llu << (height + depth));
        positions.push_back(positions.back() + cells * (height + 1));
      }
    }
    return positions;
  }

  SizeT getBitIndex(SizeT depth, SizeT height, SizeT index) const {
    return m_rowStartingPosition[depth * m_layers + height] +
           (height + 1) * (index >> (height + depth));
  }

  SizeT m_layers;
  std::vector<SizeT> m_rowStartingPosition;
  Bitvector m_bitvector;
};

class RandomCFArray {
  /**
   * The CF array for random inputs. A value v <= n of height h = lg2(v) is
   * stored in the first free cell B_d[h][i div 2^(h+d)], d = 0, ..., W.
   * Larger values, and values rejected at every depth, go to the tree T.
   * Each index keeps its height and depth (H[i], D[i]), or a bit marking that
   * its value is in T. H and D are plain packed arrays, without the Four
   * Russians encoding.
   */
public:
  explicit RandomCFArray(const std::vector<MachineWordT> &values)
      : m_size(values.size())
      , m_W(lg2(m_size - 1) + 1) // W = ceil(log n), assuming n >= 2.
      , m_paramSize(getMSb(m_W)) // Bits of a height or a depth.
      , m_B(m_size, m_W + 1)
      , m_H(m_size * m_paramSize)
      , m_D(m_size * m_paramSize)
      , m_inTree(m_size)
      , m_epochLength(std::max<SizeT>(m_size / 8, 1)) {
    for (SizeT i = 0; i < m_size; ++i) {
      insert(i, values[i]);
    }
  }

  MachineWordT get(SizeT index) const {
    if (m_inTree.getWord(index, 1) == 1) {
      return m_T.at(index);
    }
    return m_B.getWord(getParam(m_D, index), getParam(m_H, index), index);
  }

  void set(SizeT index, MachineWordT value) {
    if (m_replacements == m_epochLength) {
      rebuild();
    }
    m_replacements++;
    remove(index);
    insert(index, value);
  }

  SizeT treeSize() const { return m_T.size(); }

private:
  /**
   * Stores a value at an index that holds none.
   */
  void insert(SizeT index, MachineWordT value) {
    assert(value > 0);
    if (value <= m_size) {
      const auto height = lg2(value);
      for (SizeT depth = 0; depth <= m_W; ++depth) {
        if (m_B.getWord(depth, height, index) == 0) { // The cell is free.
          m_B.setWord(depth, height, index, value);
          setParam(m_H, index, height);
          setParam(m_D, index, depth);
          m_inTree.setWord(index, 1, 0);
          return;
        }
      }
    }
    m_T[index] = value;
    m_inTree.setWord(index, 1, 1);
  }

  /**
   * Frees the cell or the tree entry of an index. Its parameters are
   * overwritten by the next insertion.
   */
  void remove(SizeT index) {
    if (m_inTree.getWord(index, 1) == 1) {
      m_T.erase(index);
    } else {
      m_B.setWord(getParam(m_D, index), getParam(m_H, index), index, 0);
    }
  }

  /**
   * Starts a new epoch of n/8 replacements from the current values, which
   * keeps the tree small.
   */
  void rebuild() {
    std::vector<MachineWordT> values(m_size);
    for (SizeT i = 0; i < m_size; ++i) {
      values[i] = get(i);
    }
    *this = RandomCFArray(values);
  }

  SizeT getParam(const Bitvector &params, SizeT index) const {
    return params.getWord(index * m_paramSize, m_paramSize);
  }

  void setParam(Bitvector &params, SizeT index, SizeT value) {
    params.setWord(index * m_paramSize, m_paramSize, value);
  }

  SizeT m_size;
  SizeT m_W;
  SizeT m_paramSize;
  CFArrayHierarchy m_B;
  Bitvector m_H;
  Bitvector m_D;
  Bitvector m_inTree;
  Tree m_T;
  SizeT m_epochLength;
  SizeT m_replacements = 0;
};

/**
 * Draws a value whose bitlength L has Pr(L = k) = 2^-k, uniform within its
 * bitlength class.
 */
MachineWordT sampleValue(std::mt19937_64 &rng) {
  std::geometric_distribution<SizeT> heights(0.5);
  const auto height = std::min<SizeT>(heights(rng), MACHINE_WORD_BITS - 1);
  return (1llu << height) | (rng() & getLeadingBitsMask(height));
}

void test(std::vector<MachineWordT> sample, std::mt19937_64 &rng) {
  const SizeT size = sample.size();
  RandomCFArray cfarr(sample);
  for (SizeT i = 0; i < size; ++i) {
    assert(cfarr.get(i) == sample[i]);
  }
  printf("n=%llu: read test passed, tree size %llu.\n", size,
         cfarr.treeSize());

  // Fresh values at uniform indices, over several epochs.
  std::uniform_int_distribution<SizeT> indices(0, size - 1);
  for (SizeT t = 0; t < 4 * size; ++t) {
    const auto index = indices(rng);
    sample[index] = sampleValue(rng);
    cfarr.set(index, sample[index]);
  }
  for (SizeT i = 0; i < size; ++i) {
    assert(cfarr.get(i) == sample[i]);
  }
  printf("n=%llu: update test passed, tree size %llu.\n", size,
         cfarr.treeSize());
}

bool parseUInt(const char *str, uint* pOut)
{
  errno = 0;
  char *temp;
  long val = std::strtoul(str, &temp, 0);

  if (temp == str || *temp != '\0') {
    return false;
  } else if (val >= UINT_MAX) {
    printf("Overflow detected when converting %s, using UINT_MAX\n", str);
    val = UINT_MAX;
  }
  *pOut = val;
  return true;
}

int get_opt_idx(int argc, char** argv, const char* opt, const char* shopt = nullptr)
{
  for (int i = 1; i < argc; ++i) {
    if ((strcmp(opt, argv[i]) == 0) || 
        (shopt && strcmp(shopt, argv[i]) == 0)) {
      return i;
    }
  }
  return -1;
}

uint parse_options(int argc, char** argv)
{
  const int seedIdx = get_opt_idx(argc, argv, "--seed");
  const int noRandomIdx = get_opt_idx(argc, argv, "--no-random");
  if (-1 != get_opt_idx(argc, argv, "--help", "-h")) {
    printf("Usage: %s [--seed SEED | --no-random]\n", argv[0]);
    exit(2);
  } else if (-1 != noRandomIdx && -1 != seedIdx) {
    printf("Cannot use --no-random and --seed together.\n");
    exit(2);
  } else if (argc - 1 == seedIdx) {
    printf("--seed needs a number.\n");
    exit(2);
  }

  if (-1 != noRandomIdx) {
    return 0;
  }
  if (-1 != seedIdx) {
    uint result = 0;
    if (parseUInt(argv[seedIdx + 1], &result)) {
      return result;
    }
    printf("Invalid seed %s\n", argv[seedIdx + 1]);
  }
  std::random_device rd{};
  return rd();
}

int main(int argc, char** argv) {
  const auto seed = parse_options(argc, argv);
  printf("Simulation seed is: %u\n", seed);
  std::mt19937_64 rng(seed);
  for (const SizeT size : {1000, 2000, 1 << 16}) {
    std::vector<MachineWordT> sample(size);
    for (auto &value : sample) {
      value = sampleValue(rng);
    }
    test(sample, rng);
  }
  // Equal values collide at every depth, so most of them start in the tree.
  printf("Starting sanity test - multiple large values\n");
  test(std::vector<MachineWordT>(1000, 1000), rng);
  return 0;
}
