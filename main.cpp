#include "PRNG.hpp"
#include "Stribog.hpp"

int main() {
    std::string name = "ZbinyakovKiril";
    std::vector<uint8_t> seed(name.begin(), name.end());
    seed.resize(64, 0x00);

    StreebogPRNG prng(seed);

    uint64_t num1 = prng.next_int(256);
    std::cout << "Cycle 1: " << std::hex << num1 << std::endl;

    uint64_t num2 = prng.next_int(256);
    std::cout << "Cycle 2: " << std::hex << num2 << std::endl;

    uint64_t num3 = prng.next_int(256);
    std::cout << "Cycle 3: " << std::hex << num3 << std::endl;

    return 0;
}