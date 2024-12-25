import math

def Text_to_binary(Text: str) -> str:
    return  "".join(f"{ord(i):08b}" for i in Text)

def Binary_text_to_blocks(Binary_Text: str, n: int) -> list:
    Length_of_binary_text = len(Binary_Text)
    print(Length_of_binary_text)
    Check_value = 0
    if Length_of_binary_text % n == 0:
        Number_of_blocks = Length_of_binary_text // n
    else: 
        Check_value = 1
        Number_of_blocks = math.ceil(Length_of_binary_text / n)
        Blocks_with_binary_text = []
    for i in range(Number_of_blocks):
        Blocks_with_binary_text += [""]
    match Check_value:
        case 0:
            for i in range(Number_of_blocks):
                for j in range(n):
                    Blocks_with_binary_text[i] += Binary_Text[j + n * i]
        case 1:
            for i in range(Number_of_blocks - 1):
                for j in range(n):
                    Blocks_with_binary_text[i] += Binary_Text[j + n * i]
            for j in range(n - Length_of_binary_text % n - 1):
                print(i, j, n, len(Binary_Text))
                Blocks_with_binary_text[Number_of_blocks - 1] += Binary_Text[n * (i + 1) + j]
            Blocks_with_binary_text[Number_of_blocks - 1] += "0" * (n - Length_of_binary_text % n)

    return Blocks_with_binary_text

print(Text_to_binary("hello world"))

Result = Binary_text_to_blocks(Text_to_binary("hello world"), 3)
for i in range(len(Result)):
    print(Result[i])
