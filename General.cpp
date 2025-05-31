#include "Merkle_tree.hpp"
#include "PRNG.hpp"
#include "Shnorr_signature.hpp"
#include "Stribog.hpp"
#include "bigint.hpp"
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

std::vector<uint8_t>
tranformSign(std::pair<BigInt, BigInt> signature) { // R, s -> 128, 64
  std::vector<uint8_t> Rvec = signature.first.toBytes();
  std::vector<uint8_t> svec = signature.second.toBytes();
  while (Rvec.size() < 128 / 8)
    Rvec.push_back(0);
  while (svec.size() < 64 / 8)
    svec.push_back(0);

  std::vector<uint8_t> result(128 + 64, 0);
  std::copy(Rvec.begin(), Rvec.end(), result.begin());       // 128
  std::copy(svec.begin(), svec.end(), result.begin() + 128); // 64
  return result;
}

bool write_vector_to_file(const std::string &filename,
                          const std::vector<uint8_t> &data) {
  std::ofstream file(filename, std::ios::binary);
  if (!file.is_open()) {
    return false;
  }
  file.write(reinterpret_cast<const char *>(data.data()), data.size());
  return !file.fail();
}

std::vector<uint8_t> read_file_to_vector(const std::string &filename) {
  std::ifstream file(filename, std::ios::binary | std::ios::ate);
  if (!file.is_open()) {
    return {};
  }
  const std::streamsize size = file.tellg();
  file.seekg(0, std::ios::beg);
  std::vector<uint8_t> buffer(size);
  if (!file.read(reinterpret_cast<char *>(buffer.data()), size)) {
    return {};
  }

  return buffer;
}

std::vector<uint8_t> GenerateTransaction(std::string name_of_file,
                                         std::string seed) {
  std::vector<uint8_t> transaction;
  std::string name = "Sidorov Evgenii and Zbinyakov Kirill \n";
  if (name_of_file[0] == '1') {
    std::vector<uint8_t> nameBytes(name.begin(), name.end());
    std::vector<uint8_t> seedBytes(seed.begin(), seed.end());
    StreebogPRNG prng(seedBytes);
    std::vector<uint8_t> randomData = prng.next_bytes(200 - nameBytes.size());
    transaction.reserve(nameBytes.size() + randomData.size());
    transaction.insert(transaction.end(), nameBytes.begin(), nameBytes.end());
    transaction.insert(transaction.end(), randomData.begin(), randomData.end());
  } else {
    std::vector<uint8_t> seedBytes(seed.begin(), seed.end());
    StreebogPRNG prng(seedBytes);
    std::vector<uint8_t> randomData = prng.next_bytes(200);
    transaction.reserve(randomData.size());
    transaction.insert(transaction.end(), randomData.begin(), randomData.end());
  }
  write_vector_to_file(name_of_file, transaction);

  vector<uint8_t> transaction_data = read_file_to_vector(name_of_file);

  std::pair<BigInt, BigInt> keys1 = generateKeys("seed" + name_of_file);
  BigInt P1 = g.modExp(keys1.second, p);
  std::pair<BigInt, BigInt> signature =
      signMessage(q, keys1.second, P1, keys1.first, g, transaction_data);

  std::vector<uint8_t> transSign = tranformSign(signature);
  transaction.reserve(transaction.size() + transSign.size());
  transaction.insert(transaction.end(), transSign.begin(), transSign.end());
  write_vector_to_file(name_of_file + "Sign", transaction);
  std::vector<uint8_t> transaction_data2 =
      read_file_to_vector(name_of_file + "Sign");
  std::cout << "Transaction " << name_of_file << ":\n"
            << to_hex(transaction_data2) << std::endl;
  return transaction_data2;
}

int main() {
  std::vector<uint8_t> transaction1;
  std::vector<uint8_t> transaction2;
  std::vector<uint8_t> transaction3;
  std::vector<uint8_t> transaction4;
  std::vector<uint8_t> transaction5;

  transaction1 = GenerateTransaction("1.txt", "seed1");
  transaction2 = GenerateTransaction("2.txt", "seed2");
  transaction3 = GenerateTransaction("3.txt", "seed3");
  transaction4 = GenerateTransaction("4.txt", "seed4");
  transaction5 = GenerateTransaction("5.txt", "seed5");

  std::vector<std::vector<uint8_t>> dataForTree = {
      transaction1, transaction2, transaction3, transaction4, transaction5};

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
    std::cout << std::hex << std::setw(2) << std::setfill('0')
              << static_cast<int>(byte);
  }
  std::cout << std::endl;

  for (BigInt nonce = 0; nonce < 1000; nonce = nonce + 1) {
    dataForHashtry.resize(dataForHashtry.size() + 4, 0x00);
    std::vector<uint8_t> nonceBytes = nonce.toBytes();
    std::copy(nonceBytes.begin(), nonceBytes.end(),
              dataForHashtry.end() - nonceBytes.size());
    std::vector<uint8_t> hash = stribog(dataForHashtry, true);
    if ((hash[0] & 0xf8) == 0x00) {
      std::cout << "Block found! Nonce: " << nonce << std::endl;
      std::cout << "Hash: " << std::endl;
      for (unsigned char byte : hash) {
        std::cout << std::hex << std::setw(2) << std::setfill('0')
                  << static_cast<int>(byte);
      }
      std::cout << std::endl;
      break;
    }
    dataForHashtry.resize(dataForHashtry.size() - 4);
  }
}
