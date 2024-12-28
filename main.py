from GaloisMoncwari import *

p = int(input("Enter the modulus:\n"))
n = int(input("Enter the degree:\n"))
Galois = Galois_field(p, n)
poly = input("Enter the irreducible polynomial:\n")
print(Galois.Obtaining_an_irreducible_polynomial(poly))
a = input("Enter the value of first key:\n")
b = input("Enter the value of second key:\n")

mode = int(input("Enter 0 for encryption and 1 for decryption:\n"))
source = int(input("Enter 0 for plain text, 1 for binary and 2 for polynomial:\n"))
text = input("Enter the text:\n")

if mode == 0:
    if source == 0:
        print("Your text:", text)
        binary_text = Text_to_binary(text)
        print("Your binary text:", binary_text)
        blocks = Binary_text_to_blocks(binary_text, n)
        print("Your blocks:", blocks)
    elif source == 1:
        binary_text = text
        print("Your binary text:", binary_text)
        blocks = Binary_text_to_blocks(binary_text, n)
        print("Your blocks:", blocks)
    elif source == 2:
        blocks = list(text.split("/"))
    Affine_encryption(blocks, Galois, a, b)

elif mode == 1:
    if source == 0:
        print("Your text:", text)
        binary_text = Text_to_binary(text)
        print("Your binary text:", binary_text)
        blocks = Binary_text_to_blocks(binary_text, n)
        print("Your blocks:", blocks)
    elif source == 1:
        binary_text = text
        print("Your binary text:", binary_text)
        blocks = Binary_text_to_blocks(binary_text, n)
        print("Your blocks:", blocks)
    elif source == 2:
        blocks = list(text.split("/"))
    Affine_decryption(blocks, Galois, a, b)
