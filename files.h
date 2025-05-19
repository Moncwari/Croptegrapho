#pragma once

#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdint.h>
#include <vector>

void fileCleaning(std::string filename);
void writeToFile(const std::string &filename,
                 const std::vector<unsigned char> &data);

class FileByteIterator {
public:
  using BlockType = std::vector<unsigned char>;

  explicit FileByteIterator(const std::string &filename)
      : file(filename, std::ios::binary), isEnd(false) {
    if (!file) {
      throw std::runtime_error("Ошибка открытия файла");
    }
    loadNextBlock();
  }

  FileByteIterator() : isEnd(true) {}

  const BlockType &operator*() const { return block; }

  FileByteIterator &operator++() {
    if (!isEnd) {
      loadNextBlock();
    }
    return *this;
  }

  bool operator!=(const FileByteIterator &other) const { return !isEnd; }

private:
  std::ifstream file;
  BlockType block;
  bool isEnd;

  void loadNextBlock() {
    block.assign(16, 0);
    file.read(reinterpret_cast<char *>(block.data()), 16);
    std::streamsize bytesRead = file.gcount();

    if (bytesRead == 0) {
      isEnd = true;
    } else if (bytesRead < 16) {
      block.resize(bytesRead);
      isEnd = true;
    }
  }
};

class FileByteRange {
public:
  explicit FileByteRange(const std::string &filename) : filename(filename) {}

  FileByteIterator begin() const { return FileByteIterator(filename); }
  FileByteIterator end() const { return FileByteIterator(); }

private:
  std::string filename;
};

void writeToFile(const std::string &filename,
                 const std::vector<unsigned char> &data) {
  std::ofstream file(filename, std::ios::binary | std::ios::app);
  if (!file) {
    throw std::runtime_error("Ошибка открытия файла для записи");
  }
  file.write(reinterpret_cast<const char *>(data.data()), data.size());
}

void fileCleaning(std::string filename) {
  std::ofstream file(filename, std::ios::binary | std::ios::trunc);
  if (!file) {
    throw std::runtime_error("Ошибка открытия файла для очистки");
  }
}

void EncryptFromFileECB(const std::string &filenamein,
                        const std::string &filenameout,
                        const std::vector<unsigned char> &key, int len) {
  fileCleaning(filenameout);
  AES aes(AESKeyLength::AES_128);
  if (len == 192)
    aes = AES(AESKeyLength::AES_192);
  else if (len == 256)
    aes = AES(AESKeyLength::AES_256);
  long long counter = 0;
  FileByteRange fileRange(filenamein);
  std::vector<unsigned char> t;
  std::vector<unsigned char> v;
  std::vector<unsigned char> data(2048);
  for (auto it = fileRange.begin(); it != fileRange.end(); ++it) {
    v = *it;
    t = aes.EncryptECB(v, key);
    for (unsigned int i = 0; i < v.size(); i++) {
      data[counter] = t[i];
      counter++;
    }
    if (counter % 2048 == 0) {
      writeToFile(filenameout, data);
      counter = 0;
    }
  }
  std::vector<unsigned char> slice(data.begin(), data.begin() + counter);
  writeToFile(filenameout, slice);
  std::cout << "Encryption complete\n";
}

void DecryptFromFileECB(const std::string &filenamein,
                        const std::string &filenameout,
                        const std::vector<unsigned char> &key, int len) {
  fileCleaning(filenameout);
  AES aes(AESKeyLength::AES_128);
  if (len == 192)
    aes = AES(AESKeyLength::AES_192);
  else if (len == 256)
    aes = AES(AESKeyLength::AES_256);
  long long counter = 0;
  FileByteRange fileRange(filenamein);
  std::vector<unsigned char> t;
  std::vector<unsigned char> v;
  std::vector<unsigned char> data(2048);
  for (auto it = fileRange.begin(); it != fileRange.end(); ++it) {
    v = *it;
    t = aes.DecryptECB(v, key);
    for (unsigned int i = 0; i < v.size(); i++) {
      data[counter] = t[i];
      counter++;
    }
    if (counter % 2048 == 0) {
      writeToFile(filenameout, data);
      counter = 0;
    }
  }
  std::vector<unsigned char> slice(data.begin(), data.begin() + counter);
  writeToFile(filenameout, slice);
  std::cout << "Decryption complete\n";
}

void EncryptFromFileCBC(const std::string &filenamein,
                        const std::string &filenameout,
                        const std::vector<unsigned char> &key,
                        const std::vector<unsigned char> &iv, const int len) {
  fileCleaning(filenameout);
  AES aes(AESKeyLength::AES_128);
  if (len == 192)
    aes = AES(AESKeyLength::AES_192);
  else if (len == 256)
    aes = AES(AESKeyLength::AES_256);
  std::vector<uint8_t> data;
  std::vector<uint8_t> tmp;
  FileByteRange fileRange(filenamein);
  for (auto it = fileRange.begin(); it != fileRange.end(); ++it) {
    tmp = *it;
    for (unsigned int i = 0; i < tmp.size(); i++) {
      data.push_back(tmp[i]);
    }
  }

  std::vector<unsigned char> t = aes.EncryptCBC(data, key, iv);
  writeToFile(filenameout, t);
  std::cout << "Encryption complete\n";
}

void DecryptFromFileCBC(const std::string &filenamein,
                        const std::string &filenameout,
                        const std::vector<unsigned char> &key,
                        const std::vector<unsigned char> &iv, const int len) {
  fileCleaning(filenameout);
  AES aes(AESKeyLength::AES_128);
  if (len == 192)
    aes = AES(AESKeyLength::AES_192);
  else if (len == 256)
    aes = AES(AESKeyLength::AES_256);
  std::vector<uint8_t> data;
  std::vector<uint8_t> tmp;
  FileByteRange fileRange(filenamein);
  for (auto it = fileRange.begin(); it != fileRange.end(); ++it) {
    tmp = *it;
    for (unsigned int i = 0; i < tmp.size(); i++) {
      data.push_back(tmp[i]);
    }
  }

  std::vector<unsigned char> t = aes.EncryptCBC(data, key, iv);
  writeToFile(filenameout, t);
  std::cout << "Decryption complete\n";
}

void EncryptFromFileCFB(const std::string &filenamein,
                        const std::string &filenameout,
                        const std::vector<unsigned char> &key,
                        std::vector<unsigned char> iv, const int len) {
  fileCleaning(filenameout);
  AES aes(AESKeyLength::AES_128);
  if (len == 192)
    aes = AES(AESKeyLength::AES_192);
  else if (len == 256)
    aes = AES(AESKeyLength::AES_256);
  std::vector<uint8_t> data;
  std::vector<uint8_t> tmp;
  FileByteRange fileRange(filenamein);
  for (auto it = fileRange.begin(); it != fileRange.end(); ++it) {
    tmp = *it;
    for (unsigned int i = 0; i < tmp.size(); i++) {
      data.push_back(tmp[i]);
    }
  }
  std::cout << data.size() % 16 << std::endl;

  std::vector<unsigned char> t = aes.EncryptCFB(data, key, iv);
  writeToFile(filenameout, t);
  std::cout << "Encryption complete\n";
}

void DecryptFromFileCFB(const std::string &filenamein,
                        const std::string &filenameout,
                        const std::vector<unsigned char> &key,
                        std::vector<unsigned char> iv, const int len) {
  fileCleaning(filenameout);

  AES aes(AESKeyLength::AES_128);
  if (len == 192)
    aes = AES(AESKeyLength::AES_192);
  else if (len == 256)
    aes = AES(AESKeyLength::AES_256);

  std::vector<uint8_t> data;
  std::vector<uint8_t> tmp;
  FileByteRange fileRange(filenamein);

  for (auto it = fileRange.begin(); it != fileRange.end(); ++it) {
    tmp = *it;
    for (unsigned int i = 0; i < tmp.size(); i++) {
      data.push_back(tmp[i]);
    }
  }

  std::vector<unsigned char> t = aes.DecryptCFB(data, key, iv);
  writeToFile(filenameout, t);
  std::cout << "Decryption complete\n";
}
