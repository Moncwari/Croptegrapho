#include "AES.cpp"
#include "files.h"
int main() {
  std::cout << "Enter input file name: \n";
  std::string filenamein;
  std::cin >> filenamein;
  std::cout << "Enter output file name: \n";
  std::string filenameout;
  std::cin >> filenameout;
  std::cout << "Enter key length in bits: \n128\n192\n256\n";
  int key_length;
  std::cin >> key_length;
  if (key_length != 128 && key_length != 192 && key_length != 256) {
    std::cout << "Invalid key length\n";
    return -1;
  }
  std::cout << "Enter key: \n";
  std::vector<unsigned char> key;
  for (int i = 0; i < key_length / 8; i++) {
    int c;
    std::cin >> c;
    key.push_back(static_cast<unsigned char>(c));
  }
  std::cout << "Choose mode: \n1. ECB\n2. CBC\n3. CFB\n";
  int mode;
  std::cin >> mode;
  if (mode != 1 && mode != 2 && mode != 3) {
    std::cout << "Invalid mode\n";
    return -1;
  }
  std::vector<unsigned char> iv;
  if (mode != 1) {
    std::cout << "Enter IV: \n";
    for (int i = 0; i < 16; i++) {
      int c;
      std::cin >> c;
      key.push_back(static_cast<unsigned char>(c));
    }
  }

  std::cout << "Choose operation: \n1. Encryption\n2. Decryption\n";
  int operation;
  std::cin >> operation;
  if (operation == 1) {
    if (mode == 1) {
      EncryptFromFileECB(filenamein, filenameout, key, key_length);
    } else if (mode == 2) {
      EncryptFromFileCBC(filenamein, filenameout, key, iv, key_length);
    } else if (mode == 3) {
      EncryptFromFileCFB(filenamein, filenameout, key, iv, key_length);
    }
  } else if (operation == 2) {
    if (mode == 1) {
      DecryptFromFileECB(filenamein, filenameout, key, key_length);
    } else if (mode == 2) {
      DecryptFromFileCBC(filenamein, filenameout, key, iv, key_length);
    } else if (mode == 3) {
      DecryptFromFileCFB(filenamein, filenameout, key, iv, key_length);
    }
  }
  return 0;
}
