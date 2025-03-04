#include "AES.cpp"

int main() {

    std::vector<unsigned char> key = {0x2b, 0x7e, 0x15, 0x16, 0x28, 0xae, 0xd2, 0xa6, 0xab, 0xf7, 0x15, 0x88, 0x09, 0xcf, 0x4f, 0x3c};
    std::vector<unsigned char> iv = {0x2b, 0x7e, 0x15, 0x16, 0x28, 0xae, 0xd2, 0xa6, 0xab, 0xf7, 0x15, 0x88, 0x09, 0xcf, 0x4f, 0x3c};
    std::vector<unsigned char> in = {0x6b, 0xc1, 0xbe, 0xe2, 0x2e, 0x40, 0x9f, 0x96, 0xe9, 0x3d, 0x7e, 0x11, 0x73, 0x93, 0x17, 0x2a};

    AES aes(AESKeyLength::AES_128);
    std::vector<unsigned char> out1 = aes.EncryptCFB(in, key, iv);
    std::vector<unsigned char> out2 = aes.DecryptCFB(out1, key, iv);

    for (unsigned int i = 0; i < in.size(); i++) {
        if (in[i] != out2[i]) {
            std::cout << "Error\n";
        }
    }

    std::vector<unsigned char> out3 = aes.EncryptCBC(in, key, iv);
    std::vector<unsigned char> out4 = aes.DecryptCBC(out3, key, iv);

    for (unsigned int i = 0; i < in.size(); i++) {
        if (in[i] != out4[i]) {
            std::cout << "Error\n";
        }
    }

    return 0;

}