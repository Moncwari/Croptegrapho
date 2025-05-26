#pragma once

#include "Stribog.hpp"
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

std::vector<uint8_t> stribog(const std::vector<uint8_t> &message,
                             bool output_256);

class StreebogPRNG {
private:
  std::vector<uint8_t> h0_;
  uint64_t i_;

public:
  explicit StreebogPRNG(const std::vector<uint8_t> &seed) {
    if (seed.size() != 64) {
      throw std::invalid_argument("Seed must be 512 bits (64 bytes)");
    }
    h0_ = stribog(seed, true);
    i_ = 1;
  }

  std::vector<uint8_t> next_bytes(size_t length) {
    std::vector<uint8_t> output;
    while (output.size() < length) {
      std::vector<uint8_t> input_data = h0_;

      std::vector<uint8_t> i_bytes(32, 0);
      for (int j = 31; j >= 0; --j) {
        i_bytes[j] = (i_ >> (8 * (31 - j))) & 0xFF;
      }
      input_data.insert(input_data.end(), i_bytes.begin(), i_bytes.end());

      std::vector<uint8_t> hi = stribog(input_data, true);

      output.insert(output.end(), hi.begin(), hi.end());

      ++i_;
      if (i_ == 0)
        throw std::overflow_error("PRNG overflow");
    }
    output.resize(length);
    return output;
  }

  uint64_t next_int(uint32_t bits = 256) {
    if (bits == 0 || bits > 256) {
      throw std::invalid_argument("Bits must be between 1 and 256");
    }
    size_t byte_len = (bits + 7) / 8;
    auto bytes = next_bytes(byte_len);

    uint64_t result = 0;
    for (uint8_t byte : bytes) {
      result = (result << 8) | byte;
    }
    result >>= (byte_len * 8 - bits);
    return result;
  }

  void force_set_counter(uint64_t new_i) { i_ = new_i; }
};
