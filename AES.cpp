#include "AES.h"


int main() {
    vector<uint8_t> key = { 0x2b, 0x7e, 0x15, 0x16, 0x28, 0xae, 0xd2, 0xa6, 0xab, 0xf7, 0x7f, 0x67, 0x98, 0x98, 0x98, 0x98 };
    SetKeySize(128);
    vector<uint32_t> roundKeys(Nb * (Nr + 1));
    
    array<array<uint8_t, 4>, 4> state = {{{0x32, 0x43, 0xf6, 0xa8},
                                           {0x88, 0x5a, 0x30, 0x8d},
                                           {0x31, 0x31, 0x98, 0xa2},
                                           {0xe0, 0x37, 0x07, 0x34}}};
    
    cout << "Исходный блок:" << endl;
    PrintState(state);
    
    EncryptBlock(state, roundKeys);
    cout << "Зашифрованный блок:" << endl;
    PrintState(state);
    
    DecryptBlock(state, roundKeys);
    cout << "Расшифрованный блок:" << endl;
    PrintState(state);
    
    return 0;
}