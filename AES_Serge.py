from sympy.abc import x
import sympy
from itertools import product
from sympy import Poly, Pow, sympify, degree
import math
from tabulate import tabulate

S_box = [
    ["63", "7c", "77", "7b", "f2", "6b", "6f", "c5", "30", "01", "67", "2b", "fe", "d7", "ab", "76"],
    ["ca", "82", "c9", "7d", "fa", "59", "47", "f0", "ad", "d4", "a2", "af", "9c", "a4", "72", "c0"],
    ["b7", "fd", "93", "26", "36", "3f", "f7", "cc", "34", "a5", "e5", "f1", "71", "d8", "31", "15"],
    ["04", "c7", "23", "c3", "18", "96", "05", "9a", "07", "12", "80", "e2", "eb", "27", "b2", "75"],
    ["09", "83", "2c", "1a", "1b", "6e", "5a", "a0", "52", "3b", "d6", "b3", "29", "e3", "2f", "84"],
    ["53", "d1", "00", "ed", "20", "fc", "b1", "5b", "6a", "cb", "be", "39", "4a", "4c", "58", "cf"],
    ["d0", "ef", "aa", "fb", "43", "4d", "33", "85", "45", "f9", "02", "7f", "50", "3c", "9f", "a8"],
    ["51", "a3", "40", "8f", "92", "9d", "38", "f5", "bc", "b6", "da", "21", "10", "ff", "f3", "d2"],
    ["cd", "0c", "13", "ec", "5f", "97", "44", "17", "c4", "a7", "7e", "3d", "64", "5d", "19", "73"],
    ["60", "81", "4f", "dc", "22", "2a", "90", "88", "46", "ee", "b8", "14", "de", "5e", "0b", "db"],
    ["e0", "32", "3a", "0a", "49", "06", "24", "5c", "c2", "d3", "ac", "62", "91", "95", "e4", "79"],
    ["e7", "c8", "37", "6d", "8d", "d5", "4e", "a9", "6c", "56", "f4", "ea", "65", "7a", "ae", "08"],
    ["ba", "78", "25", "2e", "1c", "a6", "b4", "c6", "e8", "dd", "74", "1f", "4b", "bd", "8b", "8a"],
    ["70", "3e", "b5", "66", "48", "03", "f6", "0e", "61", "35", "57", "b9", "86", "c1", "1d", "9e"],
    ["e1", "f8", "98", "11", "69", "d9", "8e", "94", "9b", "1e", "87", "e9", "ce", "55", "28", "df"],
    ["8c", "a1", "89", "0d", "bf", "e6", "42", "68", "41", "99", "2d", "0f", "b0", "54", "bb", "16"]
    ]

Inv_S_box = [
    ["52", "09", "6a", "d5", "30", "36", "a5", "38", "bf", "40", "a3", "9e", "81", "f3", "d7", "fb"],
    ["7c", "e3", "39", "82", "9b", "2f", "ff", "87", "34", "8e", "43", "44", "c4", "de", "e9", "cb"],
    ["54", "7b", "94", "32", "a6", "c2", "23", "3d", "ee", "4c", "95", "0b", "42", "fa", "c3", "4e"],
    ["08", "2e", "a1", "66", "28", "d9", "24", "b2", "76", "5b", "a2", "49", "6d", "8b", "d1", "25"],
    ["72", "f8", "f6", "64", "86", "68", "98", "16", "d4", "a4", "5c", "cc", "5d", "65", "b6", "92"],
    ["6c", "70", "48", "50", "fd", "ed", "b9", "da", "5e", "15", "46", "57", "a7", "8d", "9d", "84"],
    ["90", "d8", "ab", "00", "8c", "bc", "d3", "0a", "f7", "e4", "58", "05", "b8", "b3", "45", "06"],
    ["d0", "2c", "1e", "8f", "ca", "3f", "0f", "02", "c1", "af", "bd", "03", "01", "13", "8a", "6b"],
    ["3a", "91", "11", "41", "4f", "67", "dc", "ea", "97", "f2", "cf", "ce", "f0", "b4", "e6", "73"],
    ["96", "ac", "74", "22", "e7", "ad", "35", "85", "e2", "f9", "37", "e8", "1c", "75", "df", "6e"],
    ["47", "f1", "1a", "71", "1d", "29", "c5", "89", "6f", "b7", "62", "0e", "aa", "18", "be", "1b"],
    ["fc", "56", "3e", "4b", "c6", "d2", "79", "20", "9a", "db", "c0", "fe", "78", "cd", "5a", "f4"],
    ["1f", "dd", "a8", "33", "88", "07", "c7", "31", "b1", "12", "10", "59", "27", "80", "ec", "5f"],
    ["60", "51", "7f", "a9", "19", "b5", "4a", "0d", "2d", "e5", "7a", "9f", "93", "c9", "9c", "ef"],
    ["a0", "e0", "3b", "4d", "ae", "2a", "f5", "b0", "c8", "eb", "bb", "3c", "83", "53", "99", "61"],
    ["17", "2b", "04", "7e", "ba", "77", "d6", "26", "e1", "69", "14", "63", "55", "21", "0c", "7d"]
    ]

Rcon = [
    ["01", "00", "00", "00"],
    ["02", "00", "00", "00"],
    ["04", "00", "00", "00"],
    ["08", "00", "00", "00"],
    ["10", "00", "00", "00"],
    ["20", "00", "00", "00"],
    ["40", "00", "00", "00"],
    ["80", "00", "00", "00"],
    ["1b", "00", "00", "00"],
    ["36", "00", "00", "00"]
    ]

class GaloisField:
    def __init__(self, p, n):  # (p-1) - макс.коэф, n - макс.степ
        self.p = p
        self.n = n

    def math_oper(self, first, second, oper, poly1):
        self.first = Poly(sympify(first), x).set_modulus(self.p)
        self.second = Poly(sympify(second), x).set_modulus(self.p)
        self.oper = oper
        self.poly1 = Poly(sympify(poly1), x).set_modulus(self.p)
        if self.oper == "+":
            ans = Poly(self.first + self.second, x).set_modulus(self.p)
            ans = Poly([coeff % self.p for coeff in ans.all_coeffs()], x)
            return ans
        if self.oper == "-":
            ans = Poly(self.first - self.second, x).set_modulus(self.p)
            ans = Poly([coeff % self.p for coeff in ans.all_coeffs()], x)
            return ans
        if self.oper == "*":
            ans = Poly(self.first * self.second, x)
            ans = Poly(ans, x).div(Poly(self.poly1, x))[1]
            ans = Poly([coeff % self.p for coeff in ans.all_coeffs()], x)
            return ans
        if self.oper == ".":
            first = Poly(self.first, x).all_coeffs()
            second = Poly(self.second, x).all_coeffs()
            while len(first) != 8:
                first.insert(0, 0)
            while len(second) != 8:
                second.insert(0, 0)
            ans = Poly.from_list([(first[i]+second[i])%2 for i in range (8)], gens=x)
            ans = Poly([coeff % self.p for coeff in ans.all_coeffs()], x)
            return ans

def to_from_block(iterable, key): # key = 0: text -> blocks; key = 1: blocks -> text
    if key == 0:
        hex_list = [hex(int(ord(x)))[2:] for x in iterable]
        while len(hex_list) != 16:
            hex_list.append("2e")
        block = [[], [], [], []]
        for i in range (0, 16):
            num = hex_list[i]
            block[i%4].append("0"*(2-len(num)) + num)
        return block

    if key == 1:
        string = ""
        block_string = ""
        for i in range (4):
            for ii in range (4):
                byte_1 = iterable[ii][i]
                block_string += chr(int(byte_1, 16))
        string += block_string
    return string
#a = to_from_block("abcdefghijklmnop", 0)
#print(a)
#print(to_from_block(a, 1))
##print()

def True_hex(n):
    n = hex(n)[2:]
    while len(n) < 2:
        n = "0" + n
    return n

def True_xor(hex1, hex2):
    bin_hex = {"0": "0000", "1": "0001", "2": "0010", "3": "0011",
               "4": "0100", "5": "0101", "6": "0110", "7": "0111",
               "8": "1000", "9": "1001", "a": "1010", "b": "1011",
               "c": "1100", "d": "1101", "e": "1110", "f": "1111"}
    bin_hex1 = ""
    bin_hex2 = ""
    for i in range (len(hex1)):
        bin_hex1 += bin_hex[hex1[i]]
        bin_hex2 += bin_hex[hex2[i]]
    bin_ans = ""
    for i in range (len(bin_hex1)):
        if bin_hex1[i] != bin_hex2[i]:
            bin_ans += "1"
        else:
            bin_ans += "0"
    ans = ""
    for i in range (0, len(bin_ans)-3, 4):
        four = bin_ans[i:i+4]
        index = list(bin_hex.values()).index(four)
        ans += list(bin_hex.keys())[index]
    return ans

def SubBytes(block):
    for i in range (4):
        for ii in range (4):
            byte = block[i][ii]
            block[i][ii] = S_box[int(byte[0], 16)][int(byte[1], 16)]
    return block

def ShiftRows(block):
    new_block = []
    for i in range (4):
        new_line_of_block = [0, 0, 0, 0]
        for ii in range(4):
            byte_1 = block[i][ii]
            new_line_of_block[(ii-i) % 4] = byte_1
        new_block.append(new_line_of_block)
    return new_block

def One_mix(column):

    field = GaloisField(2, 8)
    poly_1 = Poly.from_list([1], gens=x)
    poly_2 = Poly.from_list([1, 0], gens=x)
    poly_3 = Poly.from_list([1, 1], gens=x)
    polys = [poly_2, poly_3, poly_1, poly_1]
    field_poly = Poly.from_list([1, 0, 0, 0, 1, 1, 0, 1, 1], gens=x)

    new_column = []
    for i in range (4):
        mults = []
        for ii in range (4):
            mult = field.math_oper(polys[(ii-i)%4], column[ii], "*", field_poly)
            mults.append(mult)
        summary = field.math_oper(mults[0], mults[1], ".", field_poly)
        summary = field.math_oper(summary, mults[2], ".", field_poly)
        summary = field.math_oper(summary, mults[3], ".", field_poly)
        new_column.append(summary)
    for i in range (4):
        new_column[i] = new_column[i].all_coeffs()
    return new_column

def MixColumns(block):
    new_block = block
    for i in range (4):
        poly_old_column = []
        for ii in range (4):
            bin_byte = bin(int(block[ii][i], 16))
            poly_old_column.append(Poly.from_list([int(x) for x in bin_byte[2:]], gens=x))
        new_bin_column = One_mix(poly_old_column)
        bin_str_column = [[str(x) for x in new_bin_column[i]] for i in range (4)]
        hex_column = [True_hex(int("".join(bin_str_column[i]), 2)) for i in range (4)]
        for ii in range (4):
            new_block[ii][i] = hex_column[ii]
    return new_block

def RotWord(word):
    return [word[1], word[2], word[3], word[0]]

def SubWord(word):
    for i in range (4):
        byte = word[i]
        word[i] = S_box[int(byte[0], 16)][int(byte[1], 16)]
    return word

def XOR(matrix1, matrix2):
    matrix = matrix1
    for i in range (4):
        for ii in range (4):
            elem1 = matrix1[i][ii]
            elem2 = matrix2[i][ii]
            matrix[i][ii] = True_xor(elem1, elem2)
    return matrix

def AddRoundKey(key):
    Nk = len(key)//4
    w_i = []
    for i in range (Nk):
        w = []
        for ii in range (4):
            w.append(key[i*4+ii])
        w_i.append(w)

    for i in range (math.ceil((3*Nk+28)/Nk)):
        temp = w_i[-1]
        new_temp_str = True_xor("".join(SubWord(RotWord(temp))), "".join(Rcon[i]))
        w = "".join(w_i[Nk*i])
        new_w_str = True_xor(new_temp_str, w)
        new_w = [new_w_str[i]+new_w_str[i+1] for i in range (0, 7, 2)]
        w_i.append(new_w)

        if Nk != 8:
            for j in range (1, Nk):
                temp = "".join(w_i[-1])
                w = "".join(w_i[Nk*i+j])
                new_w_str = True_xor(temp, w)
                new_w = [new_w_str[i]+new_w_str[i+1] for i in range (0, 7, 2)]
                w_i.append(new_w)
        else:
            for j in range (1, Nk):
                if j != 4:
                    temp = "".join(w_i[-1])
                    w = "".join(w_i[Nk*i+j])
                    new_w_str = True_xor(temp, w)
                    new_w = [new_w_str[i]+new_w_str[i+1] for i in range (0, 7, 2)]
                    w_i.append(new_w)
                else:
                    temp = w_i[-1]
                    new_temp_str = "".join(SubWord(temp.copy()))
                    w = "".join(w_i[Nk*i+4])
                    new_w_str = True_xor(new_temp_str, w)
                    new_w = [new_w_str[i]+new_w_str[i+1] for i in range (0, 7, 2)]
                    w_i.append(new_w)

    keys = []
    i = 0
    while len(keys) < Nk + 7:
        transp_key = w_i[i:i+4]
        key = []
        for ii in range (4):
            line_of_key = []
            for iii in range (4):
                line_of_key.append(transp_key[iii][ii])
            key.append(line_of_key)
        keys.append(key)
        i += 4
    return keys

def InvSubBytes(block):
    for i in range (4):
        for ii in range (4):
            byte = block[i][ii]
            block[i][ii] = Inv_S_box[int(byte[0], 16)][int(byte[1], 16)]
    return block

def InvShiftRows(block):
    new_block = []
    for i in range (4):
        new_line_of_block = [0, 0, 0, 0]
        for ii in range(4):
            byte_1 = block[i][ii]
            new_line_of_block[(ii+i) % 4] = byte_1
        new_block.append(new_line_of_block)
    return new_block

def Inv_one_mix(column):
    field = GaloisField(2, 8)
    poly_9 = Poly.from_list([1, 0, 0, 1], gens=x)
    poly_b = Poly.from_list([1, 0, 1, 1], gens=x)
    poly_d = Poly.from_list([1, 1, 0, 1], gens=x)
    poly_e = Poly.from_list([1, 1, 1, 0], gens=x)
    polys = [poly_e, poly_b, poly_d, poly_9]
    field_poly = Poly.from_list([1, 0, 0, 0, 1, 1, 0, 1, 1], gens=x)

    new_column = []
    for i in range (4):
        mults = []
        for ii in range (4):
            mult = field.math_oper(polys[(ii-i)%4], column[ii], "*", field_poly)
            mults.append(mult)
        summary = field.math_oper(mults[0], mults[1], ".", field_poly)
        summary = field.math_oper(summary, mults[2], ".", field_poly)
        summary = field.math_oper(summary, mults[3], ".", field_poly)
        new_column.append(summary)
    for i in range (4):
        new_column[i] = new_column[i].all_coeffs()
    return new_column

def InvMixColumns(block):
    new_block = block
    for i in range (4):
        poly_old_column = []
        for ii in range (4):
            bin_byte = bin(int(block[ii][i], 16))
            poly_old_column.append(Poly.from_list([int(x) for x in bin_byte[2:]], gens=x))
        new_bin_column = Inv_one_mix(poly_old_column)
        bin_str_column = [[str(x) for x in new_bin_column[i]] for i in range (4)]
        hex_column = [True_hex(int("".join(bin_str_column[i]), 2)) for i in range (4)]
        for ii in range (4):
            new_block[ii][i] = hex_column[ii]
    return new_block

def transp(block):
    new_block = []
    for i in range (4):
        new_line = []
        for ii in range (4):
            new_line.append(block[ii][i])
        new_block.append(new_line)
    return new_block

def AES_cipher(message, key, flag):
    block_of_message = to_from_block(message, 0)
    key = key.split(" ")
    keys = AddRoundKey(key)
    Nk = len(key)//4

    if flag == 0:
        print(block_of_message)
        print(keys[0])
        block_of_message = XOR(block_of_message, keys[0])
        print("Before round 1")
        print(tabulate(block_of_message, tablefmt="grid"))
        for i in range (1, Nk+6):
            block_of_message = XOR(MixColumns(ShiftRows(SubBytes(block_of_message))), keys[i])
            print("Before round "+str(i+1))
            print(tabulate(block_of_message, tablefmt="grid"))
        block_of_message = (XOR(ShiftRows(SubBytes(block_of_message)), keys[-1]))
        print("Result")
        print(tabulate(block_of_message, tablefmt="grid"))
        return to_from_block(block_of_message, 1)
    else:
        keys = keys[::-1]
        print(block_of_message)
        print(keys[0])
        block_of_message = InvSubBytes(InvShiftRows(XOR(block_of_message, keys[0])))
        print("After round 1")
        print(tabulate(block_of_message, tablefmt="grid"))
        for i in range (1, Nk+6):
            block_of_message = InvSubBytes(InvShiftRows(InvMixColumns(XOR(block_of_message, keys[i]))))
            print("After round "+str(i+1))
            print(tabulate(block_of_message, tablefmt="grid"))
        block_of_message = XOR(block_of_message, keys[-1])
        print("Result")
        print(tabulate(block_of_message, tablefmt="grid"))
        return to_from_block(block_of_message, 1)

key = "2b 7e 15 16 28 ae d2 a6 ab f7 15 88 09 cf 4f 3c"
message = "abcdefghijklmnop"
b = AES_cipher(message, key, 0)
print()
print(b)
print()
print(AES_cipher(b, key, 1))
