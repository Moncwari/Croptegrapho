import math
import numpy as np
import random

N = 16
K = 8
L = 4

data = np.array([1, 0, 1, 0, 1, 0, 1, 0])
Tree_dict = {}
Safe_bits = ["LLLL", "LLLR", "LLRL", "LLRR", "LRLL", "LRLR"]

#etalon = [0, 0, 0, 0, 0, 0] + data[0:2] + [0, 0] + data[2:]
def concatinate(x : list, y : list):
    tmp = [x[i] ^ y[i] for i in range(len(x))]
    for i in y:
        tmp.append(i)
    return tmp

def L_func(x, y): return (x * y) // abs(x * y) * min(abs(x), abs(y))

def R_func(x, y, b):  return x + y if b == 0 else y - x

def gaussian(m : list):
    return [random.gauss(0, 2) * i for i in range(len(m))]

def BPSK(m : list):
    a = [1 if x == 0 else 0 for x in m]
    return a

def error(m : list, number=2):
    pass

def Tree_level(Elements: list, level):
    if len(Elements) == 1:
        Tree_dict[level] = Elements
        return [Elements]
    length = len(Elements)
    Result = []
    for i in range(0, length):
        #print(i)
        Result.append([Elements[i], Elements[i+length//2]])
        if i == length//2 - 1:
            break
    Tree_dict[level] = Result
    return Result

def Tree_level_L(Elements: list[list], level):
    # if len(Elements) == 1:
    #     return L_func(Elements[0][0], Elements[0][1])
    length = len(Elements)
    Help_Result = []
    for i in range(0, length):
        # print(Elements)
        Help_Result.append(L_func(Elements[i][0], Elements[i][1]))
    #print(Help_Result)
    length = len(Help_Result)
    # if length == 1:
    #     return [Help_Result]
    return Tree_level(Help_Result, level)

def Tree_level_R(Elements: list[list], level, Byte_value: int):
    length = len(Elements)
    Help_Result = []
    for i in range(0, length):
        # print(Elements)
        Help_Result.append(R_func(Elements[i][0], Elements[i][1], Byte_value))
    length = len(Help_Result)
    return Tree_level(Help_Result, level)

# def Tree():
#     level = "L"
#     Tree_dict["All"] = Elements
#     Result = Tree_level(Elements, level)
#     #print(Result)
#     length = len(Result)
#     while length >= 1:
#         Result = Tree_level_L(Result, level)
#         level += "L"
#         length = len(Result[0])  
#         if length == 1:
#             break
#     return Result

def Tree_L(Elements: list[list]):
    level = "L"
    # Tree_dict["All"] = Elements
    # Result = Tree_level(Elements, level)
    #print(Result)
    # length = len(Result)
    Result = Elements[:]
    length = len(Result)
    while length >= 1:
        Result = Tree_level_L(Result, level)
        level += "L"
        length = len(Result[0])  
        if length == 1:
            break
    return Result

# def find_leftmost_path(paths: list[str]) -> str:
#     return min(paths, key=lambda s: (s.count('R'), s))  # Сначала минимум 'R', потом лексикографически

def find_leftmost_path(paths: list[str]) -> str:
    return max(paths, key=lambda s: (s.count('L'), -s.count('R'), -len(s), ''.join(reversed(s))))


def Tree_R(Elements: list[list]):
    Help_Tree_dict = Tree_dict.copy()
    Levels = list(Help_Tree_dict.keys())
    Most_left_key = find_leftmost_path(Levels)
    value = Tree_dict[Most_left_key]
    if Most_left_key in Safe_bits:
        Byte_value = 0
    else: 
        print(Tree_dict[Most_left_key])
        if Tree_dict[Most_left_key][0][0] >= 0:
            Byte_value = 0
        else:
            Byte_value = 1

    Actual_level = Most_left_key[:-1] + "R"
    #print(list(Help_Tree_dict.keys()))
    Keys = list(Help_Tree_dict.keys())
    while len(list(Help_Tree_dict.keys())) > 1:
        Work_level = Most_left_key[:-1]
        Result = Tree_level_R(Tree_dict[Work_level], Actual_level, Byte_value)

    return Tree_dict

# 0 1 2 3 4 5 6 7 || 8 9 10 11 12 13 14 15
Example_data = [0.6, 0.7, 0.79, 0.54, 0.4, 0.8, -0.3, 0.39, 0.44, 0.36, -0.9, -0.13, -0.81, 0.62, -0.48, -0.6]
Tree_dict["all"] = Example_data
Example_data = Tree_level(Example_data, "L")


print(Tree_L(Example_data))
print(Tree_R(Example_data))
print(Tree_dict)
# print(find_leftmost_path(Tree_dict.keys()))