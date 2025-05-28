#include <vector>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <cstdint>
#include "Stribog.hpp"

using namespace std;

class MerkleTree {
private:
    vector<vector<vector<uint8_t>>> levels;
    
    vector<uint8_t> combineHashes(const vector<uint8_t>& left, const vector<uint8_t>& right) {
        vector<uint8_t> combined;
        combined.reserve(left.size() + right.size());
        combined.insert(combined.end(), left.begin(), left.end());
        combined.insert(combined.end(), right.begin(), right.end());
        return stribog(combined, true);
    }

public:
    explicit MerkleTree(const vector<vector<uint8_t>>& data) {
        if (data.empty()) {
            throw invalid_argument("Input data cannot be empty");
        }

        vector<vector<uint8_t>> leaves;
        for (const auto& block : data) {
            leaves.push_back(stribog(block, true));
        }

        vector<vector<uint8_t>> currentLevel = leaves;
        levels.push_back(currentLevel);

        while (currentLevel.size() > 1) {
            vector<vector<uint8_t>> nextLevel;
            
            for (size_t i = 0; i < currentLevel.size(); i += 2) {
                if (i + 1 >= currentLevel.size()) {
                    nextLevel.push_back(combineHashes(currentLevel[i], currentLevel[i]));
                } else {
                    nextLevel.push_back(combineHashes(currentLevel[i], currentLevel[i+1]));
                }
            }
            
            currentLevel = nextLevel;
            levels.push_back(currentLevel);
        }
    }

    vector<uint8_t> getRoot() const {
        return levels.back().front();
    }

    const vector<vector<vector<uint8_t>>>& getLevels() const {
        return levels;
    }
};