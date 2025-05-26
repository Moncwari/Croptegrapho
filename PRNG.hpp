#pragma once

#include "Stribog.hpp"
#include <cstdint>
#include <vector>
#include <stdexcept>

std::vector<uint8_t> stribog(const std::vector<uint8_t> &message, bool output_256);

class StreebogPRNG {
private:
    std::vector<uint8_t> h0_;
    uint64_t counter_;
    std::vector<uint8_t> buffer_;
    size_t buffer_offset_;

    void refill_buffer() {
        std::vector<uint8_t> input = h0_;

        std::vector<uint8_t> counter_bytes(32, 0);
        for (int j = 0; j < 32; ++j) {
            counter_bytes[31 - j] = (counter_ >> (8 * j)) & 0xFF;
        }
        input.insert(input.end(), counter_bytes.begin(), counter_bytes.end());

        buffer_ = stribog(input, true);
        buffer_offset_ = 0;

        if (++counter_ == 0) {
            throw std::overflow_error("PRNG overflow");
        }
    }

public:
    explicit StreebogPRNG(const std::vector<uint8_t> &seed) : counter_(1), buffer_offset_(0) {
        if (seed.size() != 64) {
            throw std::invalid_argument("Seed must be 512 bits (64 bytes)");
        }
        h0_ = stribog(seed, true);
        refill_buffer();
    }

    uint8_t next_byte() {
        if (buffer_offset_ >= buffer_.size()) {
            refill_buffer();
        }
        return buffer_[buffer_offset_++];
    }

    std::vector<uint8_t> next_bytes(size_t length) {
        std::vector<uint8_t> output;
        output.reserve(length);
        while (output.size() < length) {
            output.push_back(next_byte());
        }
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
};
