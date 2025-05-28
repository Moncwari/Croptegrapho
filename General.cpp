#include <vector>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <cstdint>
#include "Stribog.hpp"
#include "PRNG.hpp"
#include "bigint.hpp"
#include "Shnorr_signature.hpp"
#include "Merkle_tree.hpp"

std::vector<uint8_t> tranformSign(std::pair<BigInt, BigInt> signature) { // R, s -> 128, 64
    std::vector<uint8_t> Rvec = signature.first.toBytes();
    std::vector<uint8_t> svec = signature.second.toBytes();

    while (Rvec.size() < 128) Rvec.push_back(0);
    while (svec.size() < 64) svec.push_back(0);

    std::vector<uint8_t> result(128 + 64, 0);
    std::copy(Rvec.begin(), Rvec.end(), result.begin()); // 128
    std::copy(svec.begin(), svec.end(), result.begin() + 128); // 64

    return result;
}

int main() {
    std::vector<uint8_t> transaction1;
    std::vector<uint8_t> transaction2;
    std::vector<uint8_t> transaction3;
    std::vector<uint8_t> transaction4;
    std::vector<uint8_t> transaction5;

    std::string name = "KirillZbinyakov+AnastasiaTuturova";
    std::vector<uint8_t> nameBytes(name.begin(), name.end());

    std::string seed = "MakeLoveNotCrypto";
    std::vector<uint8_t> seedBytes(seed.begin(), seed.end());

    StreebogPRNG prng(seedBytes);
    std::vector<uint8_t> randomData1 = prng.next_bytes(200 - nameBytes.size());
    std::vector<uint8_t> randomData2 = prng.next_bytes(200);
    std::vector<uint8_t> randomData3 = prng.next_bytes(200);
    std::vector<uint8_t> randomData4 = prng.next_bytes(200);
    std::vector<uint8_t> randomData5 = prng.next_bytes(200);

    transaction1.insert(transaction1.begin(), nameBytes.begin(), nameBytes.end());
    transaction1.insert(transaction1.end(), randomData1.begin(), randomData1.end());

    transaction2.insert(transaction2.begin(), randomData2.begin(), randomData2.end());
    transaction3.insert(transaction2.begin(), randomData2.begin(), randomData2.end());
    transaction4.insert(transaction2.begin(), randomData2.begin(), randomData2.end());
    transaction5.insert(transaction2.begin(), randomData2.begin(), randomData2.end());

    std::pair<BigInt, BigInt> keys1 = generateKeys("seed1"); // (x, r)
    std::pair<BigInt, BigInt> keys2 = generateKeys("seed2");
    std::pair<BigInt, BigInt> keys3 = generateKeys("seed3");
    std::pair<BigInt, BigInt> keys4 = generateKeys("seed4");
    std::pair<BigInt, BigInt> keys5 = generateKeys("seed5");

    BigInt P1 = g ^ keys1.first % p;
    BigInt P2 = g ^ keys2.first % p;
    BigInt P3 = g ^ keys3.first % p;
    BigInt P4 = g ^ keys4.first % p;
    BigInt P5 = g ^ keys5.first % p;

    std::pair<BigInt, BigInt> signature1 = signMessage(q, keys1.second, P1, keys1.first, g, transaction1);
    std::pair<BigInt, BigInt> signature2 = signMessage(q, keys2.second, P2, keys2.first, g, transaction2);
    std::pair<BigInt, BigInt> signature3 = signMessage(q, keys3.second, P3, keys3.first, g, transaction3);
    std::pair<BigInt, BigInt> signature4 = signMessage(q, keys4.second, P4, keys4.first, g, transaction4);
    std::pair<BigInt, BigInt> signature5 = signMessage(q, keys5.second, P5, keys5.first, g, transaction5);

    transaction1.insert(transaction1.end(), tranformSign(signature1).begin(), tranformSign(signature1).end());
    transaction2.insert(transaction2.end(), tranformSign(signature2).begin(), tranformSign(signature2).end());
    transaction3.insert(transaction3.end(), tranformSign(signature3).begin(), tranformSign(signature3).end());
    transaction4.insert(transaction4.end(), tranformSign(signature4).begin(), tranformSign(signature4).end());
    transaction5.insert(transaction5.end(), tranformSign(signature5).begin(), tranformSign(signature5).end());

    std::vector<std::vector<uint8_t>> dataForTree = {transaction1, transaction2, transaction3, transaction4, transaction5};

    MerkleTree tree(dataForTree);
    std::vector<uint8_t> root = tree.getRoot();

    std::string blockSeed = "TestovaOneLove";
    std::vector<uint8_t> blockSeedBytes(blockSeed.begin(), blockSeed.end());

    StreebogPRNG prng2(blockSeedBytes);
    std::vector<uint8_t> prevHash = prng2.next_bytes(32);

    std::vector<uint8_t> timestamp = {17, 28, 05, 25};
    std::vector<uint8_t> dataForHash = {};
    dataForHash.insert(dataForHash.end(), prevHash.begin(), prevHash.end());
    dataForHash.insert(dataForHash.end(), root.begin(), root.end());
    dataForHash.insert(dataForHash.end(), timestamp.begin(), timestamp.end());
    
    for (BigInt nonce = 0; nonce < 1000000000; nonce = nonce + 1) {
        if (nonce % 1 == 0) std::cout << "Nonce: " << nonce << std::endl;
        std::vector<uint8_t> dataForHashtry = dataForHash;
        dataForHashtry.insert(dataForHashtry.end(), nonce.toBytes().begin(), nonce.toBytes().end());
        std::vector<uint8_t> hash = stribog(dataForHash, true);
        if (hash[0] & 0xf8 == 0x00) {
            std::cout << "Block found! Nonce: " << nonce << std::endl;
            std::cout << "Hash: " << std::hex << std::setfill('0') << std::setw(2 * hash.size()) << std::string(hash.begin(), hash.end()) << std::endl;
        }
    }

}
