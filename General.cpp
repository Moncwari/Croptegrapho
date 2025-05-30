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
    while (Rvec.size() < 128 / 8) Rvec.push_back(0);
    while (svec.size() < 64 / 8) svec.push_back(0);

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

    std::cout << "Random data completed\n";

    transaction1.reserve(nameBytes.size() + randomData1.size());
    transaction2.reserve(randomData2.size());
    transaction3.reserve(randomData3.size());
    transaction4.reserve(randomData4.size());
    transaction5.reserve(randomData5.size());

    transaction1.insert(transaction1.end(), nameBytes.begin(), nameBytes.end());
    transaction1.insert(transaction1.end(), randomData1.begin(), randomData1.end());

    transaction2.insert(transaction2.end(), randomData2.begin(), randomData2.end());
    transaction3.insert(transaction3.end(), randomData3.begin(), randomData3.end());
    transaction4.insert(transaction4.end(), randomData4.begin(), randomData4.end());
    transaction5.insert(transaction5.end(), randomData5.begin(), randomData5.end());
    
    std::cout << "Transactions completed\n";

    std::pair<BigInt, BigInt> keys1 = generateKeys("seed1"); // (x, r)
    std::pair<BigInt, BigInt> keys2 = generateKeys("seed2");
    std::pair<BigInt, BigInt> keys3 = generateKeys("seed3");
    std::pair<BigInt, BigInt> keys4 = generateKeys("seed4");
    std::pair<BigInt, BigInt> keys5 = generateKeys("seed5");

    std::cout << "Keys completed\n";

    BigInt P1 = g.modExp(keys1.second, p);
    BigInt P2 = g.modExp(keys2.second, p);
    BigInt P3 = g.modExp(keys3.second, p);
    BigInt P4 = g.modExp(keys4.second, p);
    BigInt P5 = g.modExp(keys5.second, p);

    std::cout << "Public keys completed\n";

    std::pair<BigInt, BigInt> signature1 = signMessage(q, keys1.second, P1, keys1.first, g, transaction1);
    std::pair<BigInt, BigInt> signature2 = signMessage(q, keys2.second, P2, keys2.first, g, transaction2);
    std::pair<BigInt, BigInt> signature3 = signMessage(q, keys3.second, P3, keys3.first, g, transaction3);
    std::pair<BigInt, BigInt> signature4 = signMessage(q, keys4.second, P4, keys4.first, g, transaction4);
    std::pair<BigInt, BigInt> signature5 = signMessage(q, keys5.second, P5, keys5.first, g, transaction5);
    
    std::cout << "Signatures completed\n";
    std::vector<uint8_t> transSign = tranformSign(signature1);
    transaction1.reserve(transaction1.size() + transSign.size());
    transaction1.insert(transaction1.end(), transSign.begin(), transSign.end());
    transSign = tranformSign(signature2);
    transaction2.reserve(transaction2.size() + transSign.size());
    transaction2.insert(transaction2.end(), transSign.begin(), transSign.end());
    transSign = tranformSign(signature3);
    transaction3.reserve(transaction3.size() + transSign.size());
    transaction3.insert(transaction3.end(), transSign.begin(), transSign.end());
    transSign = tranformSign(signature4);
    transaction4.reserve(transaction4.size() + transSign.size());
    transaction4.insert(transaction4.end(), transSign.begin(), transSign.end());
    transSign = tranformSign(signature5);
    transaction5.reserve(transaction5.size() + transSign.size());
    transaction5.insert(transaction5.end(), transSign.begin(), transSign.end());

    std::cout << "Transactions completed\n";

    std::vector<std::vector<uint8_t>> dataForTree = {transaction1, transaction2, transaction3, transaction4, transaction5};

    MerkleTree tree(dataForTree);
    std::vector<uint8_t> root = tree.getRoot();

    std::cout << "Root: \n" << to_hex(root) << std::endl;

    std::string blockSeed = "TestovaOneLove";
    std::vector<uint8_t> blockSeedBytes(blockSeed.begin(), blockSeed.end());

    StreebogPRNG prng2(blockSeedBytes);
    std::vector<uint8_t> prevHash = prng2.next_bytes(32);

    std::vector<uint8_t> timestamp = {17, 28, 05, 25};
    std::vector<uint8_t> dataForHash = {};
    dataForHash.insert(dataForHash.end(), prevHash.begin(), prevHash.end());
    dataForHash.insert(dataForHash.end(), root.begin(), root.end());
    dataForHash.insert(dataForHash.end(), timestamp.begin(), timestamp.end());
    
    std::cout << "Data for hash completed\n";
    std::vector<uint8_t> dataForHashtry = dataForHash;

    for (unsigned char byte : dataForHash) {
        std::cout << std::hex 
                << std::setw(2) 
                << std::setfill('0') 
                << static_cast<int>(byte);
    }
    std::cout << std::endl;
    
    
    for (BigInt nonce = 0; nonce < 1000; nonce = nonce + 1) {
        dataForHashtry.resize(dataForHashtry.size() + 4, 0x00);
        std::vector<uint8_t> nonceBytes = nonce.toBytes();
        std::copy(nonceBytes.begin(), nonceBytes.end(), dataForHashtry.end() - nonceBytes.size());
        std::vector<uint8_t> hash = stribog(dataForHashtry, true);
        if ((hash[0] & 0xf8) == 0x00) {
            std::cout << "Block found! Nonce: " << nonce << std::endl;
            std::cout << "Hash: " << std::endl;
            for (unsigned char byte : hash) {
                std::cout << std::hex 
                        << std::setw(2) 
                        << std::setfill('0') 
                        << static_cast<int>(byte);
            }
            std::cout << std::endl;
            break;
        }
        dataForHashtry.resize(dataForHashtry.size() - 4);
    }

}
