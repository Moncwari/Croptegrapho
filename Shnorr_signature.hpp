#pragma once

#include <iostream>
#include <string>
#include <sstream>
#include <random>
#include <bitset>
#include "bigint.hpp"
#include "PRNG.hpp"


static std::bitset<1024> hexToBitset(const std::string& hex) {
    std::bitset<1024> bs;
    int bitpos = 0;

    for (int i = (int)hex.size() - 1; i >= 0; --i) {
        char c = hex[i];
        int v = 0;
        if      (c >= '0' && c <= '9') v = c - '0';
        else if (c >= 'A' && c <= 'F') v = c - 'A' + 10;
        else if (c >= 'a' && c <= 'f') v = c - 'a' + 10;
        else continue;

        for (int b = 0; b < 4; ++b) {
            bs[bitpos++] = (v >> b) & 1;
        }
    }
    return bs;
}

const BigInt p = BigInt::bitsToBigInt(
    hexToBitset("EE8172AE8996608FB69359B89EB82A69854510E2977A4D63BC97322CE5DC3386EA0A12B343E9190F23177539845839786BB0C345D165976EF2195EC9B1C379E3")
);
const BigInt q = BigInt::bitsToBigInt(
    hexToBitset("98915E7EC8265EDFCDA31E88F24809DDB064BDC7285DD50D7289F0AC6F49DD2D")
);
const BigInt g =  BigInt::bitsToBigInt(
    hexToBitset("9E96031500C8774A869582D4AFDE2127AFAD2538B4B6270A6F7C8837B50D50F206755984A49E509304D648BE2AB5AAB18EBE2CD46AC3D8495B142AA6CE23E21C")
);

BigInt generateRandomGroupNumber(BigInt q, std::string seed = "Lebedeva") {
    StreebogPRNG prng(std::vector<uint8_t>(seed.begin(), seed.end()));
    std::vector<uint8_t> init = prng.next_bytes(63);

    StreebogPRNG prng2(init);
    std::vector<uint8_t> result = prng2.next_bytes(63);
    return BigInt::BytesToBigInt(result);
}

std::pair<BigInt, BigInt> generateKeys(std::string seed = "Lebedeva") {
    BigInt x = generateRandomGroupNumber(q);
    BigInt r = generateRandomGroupNumber(q);
    return std::make_pair(x, r);
}

std::pair<BigInt, BigInt> signMessage(BigInt q, BigInt r, BigInt P, BigInt x, BigInt g, std::vector<uint8_t> message) {
    BigInt R = g.modExp(r, p);
    std::vector<uint8_t> Rvec = R.toBytes();
    std::vector<uint8_t> Pvec = P.toBytes();
    std::vector<uint8_t> eBase;
    eBase.reserve(Rvec.size() + Pvec.size() + message.size());
    
    eBase.insert(eBase.end(), Rvec.begin(), Rvec.end());
    eBase.insert(eBase.end(), Pvec.begin(), Pvec.end());
    eBase.insert(eBase.end(), message.begin(), message.end());
    BigInt e = BigInt::BytesToBigInt(eBase);

    BigInt s = (r + e * x) % q;

    return std::make_pair(R, s);
}

bool checkSign(BigInt R, BigInt s, BigInt P, std::vector<uint8_t> message) {
    std::vector<uint8_t> Rvec = R.toBytes();
    std::vector<uint8_t> Pvec = P.toBytes();
    std::vector<uint8_t> eBase;
    eBase.reserve(Rvec.size() + Pvec.size() + message.size());
    
    eBase.insert(eBase.end(), Rvec.begin(), Rvec.end());
    eBase.insert(eBase.end(), Pvec.begin(), Pvec.end());
    eBase.insert(eBase.end(), message.begin(), message.end());
    
    BigInt e = BigInt::BytesToBigInt(eBase);

    BigInt left = R * (P^e);
    BigInt right = g^s;

    return left == right;
}