#pragma once

#include "Stribog.hpp"
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace std;

class MerkleTree {
private:
  vector<vector<vector<uint8_t>>> levels;

  vector<uint8_t> combineHashes(vector<uint8_t> left, vector<uint8_t> right) {
    vector<uint8_t> combined = std::move(left);
    combined.reserve(combined.size() + right.size());
    combined.insert(combined.end(), std::make_move_iterator(right.begin()),
                    std::make_move_iterator(right.end()));

    return stribog(combined, true);
  }

public:
  explicit MerkleTree(const vector<vector<uint8_t>> &data) {
    if (data.empty())
      throw invalid_argument("Input data cannot be empty");

    vector<vector<uint8_t>> leaves;
    leaves.reserve(data.size());

    for (const auto &block : data) {
      leaves.push_back(stribog(block, true));
    }

    levels.push_back(leaves);
    vector<vector<uint8_t>> currentLevel = std::move(leaves);

    while (currentLevel.size() > 1) {
      vector<vector<uint8_t>> nextLevel;
      nextLevel.reserve((currentLevel.size() + 1) / 2);

      for (size_t i = 0; i < currentLevel.size();) {
        if (i + 1 >= currentLevel.size()) {
          vector<uint8_t> copy = currentLevel[i];
          nextLevel.push_back(
              combineHashes(std::move(currentLevel[i]), std::move(copy)));
          i += 1;
        } else {
          nextLevel.push_back(combineHashes(std::move(currentLevel[i]),
                                            std::move(currentLevel[i + 1])));
          i += 2;
        }
      }

      levels.push_back(nextLevel);
      currentLevel = std::move(nextLevel);
    }
  }

  vector<uint8_t> getRoot() const { return levels.back().front(); }

  const vector<vector<vector<uint8_t>>> &getLevels() const { return levels; }
};