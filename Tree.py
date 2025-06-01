import math
import numpy as np
import random
import copy
import collections

N = 16
K = 8
L = 4

SAFE_BITS = ["LLLL", "LLLR", "LLRL", "LLRR", "LRLL", "LRLR", "RLLL", "RLLR"]
PATH_list = [
    "L",
    "LL",
    "LLL",
    "LLLL",
    "LLLR",
    "LLL",
    "LLR",
    "LLRL",
    "LLRR",
    "LLR",
    "LL",
    "LR",
    "LRL",
    "LRLL",
    "LRLR",
    "LRL",
    "LRR",
    "LRRL",
    "LRRR",
    "LRR",
    "LR",
    "L",
    "R",
    "RL",
    "RLL",
    "RLLL",
    "RLLR",
    "RLL",
    "RLR",
    "RLRL",
    "RLRR",
    "RLR",
    "RL",
    "RR",
    "RRL",
    "RRLL",
    "RRLR",
    "RRL",
    "RRR",
    "RRRL",
    "RRRR",
    "RRR",
    "RR",
    "R",
    "Finish",
]


def concatinate(x: list, y: list):
    tmp = [x[i] ^ y[i] for i in range(len(x))]
    for i in y:
        tmp.append(i)
    return tmp


def sign(x):
    return 0 if x >= 0 else 1


def L_func(x, y):
    return (1 if x * y >= 0 else -1) * min(abs(x), abs(y))


def R_func(x, y, b):
    return x + y if b == 0 else y - x


def gaussian(m: list):
    return [random.gauss(0, 2) * i for i in range(len(m))]


def BPSK(m: list):
    a = [1 if x == 0 else -1 for x in m]
    return a

def add_lists_elementwise(list1, list2):
    if len(list1) != len(list2):
        print("Ошибка: Списки должны быть одинаковой длины для поэлементного сложения.")
        return None

    result_list = []
    for i in range(len(list1)):
        result_list.append(list1[i] ^ list2[i])

    return result_list


def from_8_to_16(inp):
    u = [0] * 6 + inp[:2] + [0, 0] + inp[2:]
    k = int(math.sqrt(len(u)))
    l = [[] for _ in range(k)]
    for k1 in range(k - 1, -1, -1):
        if k1 == k - 1:
            for i in range(0, len(u), 2):
                l[k1].append([int(u[i]) ^ int(u[i + 1]), int(u[i + 1])])
        else:
            for j in range(0, len(l[k1 + 1]), 2):
                l[k1].append(
                    add_lists_elementwise(l[k1 + 1][j], l[k1 + 1][j + 1])
                    + l[k1 + 1][j + 1]
                )
    # print(f"Our list at k = {k1}: ", *l, sep = "\n")
    return l[0][0]


class Tree:
    def __init__(self, Elements: list[float]):
        self.dict_nodes = {}
        self.dict_nodes["ALL"] = (Elements, [])
        self.nodes = collections.deque(PATH_list)
        self.proceed = set()
        self.proceed.add("ALL")

    def Tree_level_L(self):
        if len(self.cur_node) == 4:
            if self.cur_node in SAFE_BITS:
                self.dict_nodes[self.cur_node] = {
                    "pred": L_func(
                        self.dict_nodes[self.cur_node[:-1]][0][0],
                        self.dict_nodes[self.cur_node[:-1]][0][1],
                    ),
                    "should": [0],
                }
                self.proceed.add(self.cur_node)
            else:
                predicted = L_func(
                    self.dict_nodes[self.cur_node[:-1]][0][0],
                    self.dict_nodes[self.cur_node[:-1]][0][1],
                )
                leaf1 = {"pred": predicted, "should": [0]}
                leaf2 = {"pred": predicted, "should": [1]}
                self.dict_nodes[self.cur_node] = leaf1
                self.proceed.add(self.cur_node)
                copied = copy.deepcopy(self)
                copied.dict_nodes[self.cur_node] = leaf2
                all_trees.append(copied)

        else:
            if self.cur_node not in self.dict_nodes.keys():
                node = self.cur_node[:-1]
                if node == "":
                    node = "ALL"
                tmp = self.dict_nodes[node][0]
                length = len(tmp)
                Help_Result = []
                for i in range(0, length // 2):
                    Help_Result.append(L_func(tmp[i], tmp[i + length // 2]))
                self.dict_nodes[self.cur_node] = (Help_Result, [])
            else:
                self.dict_nodes[self.cur_node] = (
                    self.dict_nodes[self.cur_node][0],
                    concatinate(
                        self.dict_nodes[self.cur_node + "L"][
                            1 if len(self.cur_node) != 3 else "should"
                        ],
                        self.dict_nodes[self.cur_node + "R"][
                            1 if len(self.cur_node) != 3 else "should"
                        ],
                    ),
                )

        return self.cur_node

    def Tree_level_R(self):
        if len(self.cur_node) == 4:
            node = self.cur_node[:-1] + "L"
            b = self.dict_nodes[node]["should"][0]
            if self.cur_node in SAFE_BITS:
                self.dict_nodes[self.cur_node] = {
                    "pred": R_func(
                        self.dict_nodes[self.cur_node[:-1]][0][0],
                        self.dict_nodes[self.cur_node[:-1]][0][1],
                        b,
                    ),
                    "should": [0],
                }
                self.proceed.add(self.cur_node)
                node = node[:-1]
                self.dict_nodes[node][1].append(
                    self.dict_nodes[self.cur_node[:-1] + "L"]["should"]
                )
                self.dict_nodes[node][1].append(
                    self.dict_nodes[self.cur_node[:-1] + "R"]["should"]
                )
                # self.cur_node = node
                return self.cur_node
            else:
                predicted = R_func(
                    self.dict_nodes[self.cur_node[:-1]][0][0],
                    self.dict_nodes[self.cur_node[:-1]][0][1],
                    b,
                )
                node = node[:-1]
                leaf1 = {"pred": predicted, "should": [0]}
                leaf2 = {"pred": predicted, "should": [1]}
                self.dict_nodes[self.cur_node] = leaf1
                self.proceed.add(self.cur_node)
                copied = copy.deepcopy(self)
                copied.dict_nodes[self.cur_node] = leaf2
                all_trees.append(copied)
                return self.cur_node

        else:
            if self.cur_node not in self.dict_nodes.keys():
                node = self.cur_node[:-1]
                if node == "":
                    node = "ALL"
                tmp = self.dict_nodes[node][0]
                length = len(tmp)
                Help_Result = []
                for i in range(0, length // 2):
                    Help_Result.append(
                        R_func(
                            tmp[i],
                            tmp[i + length // 2],
                            self.dict_nodes[node + "L" if node != "ALL" else "L"][1][i],
                        )
                    )
                self.dict_nodes[self.cur_node] = (Help_Result, [])
                return self.cur_node
            else:
                self.dict_nodes[self.cur_node] = (
                    self.dict_nodes[self.cur_node][0],
                    concatinate(
                        self.dict_nodes[self.cur_node + "L"][
                            1 if len(self.cur_node) != 3 else "should"
                        ],
                        self.dict_nodes[self.cur_node + "R"][
                            1 if len(self.cur_node) != 3 else "should"
                        ],
                    ),
                )
                return self.cur_node


def calculate_metrics(Tree: Tree):
    decode = ""
    Error = 0
    All_nodes = list(Tree.dict_nodes.keys())
    for node in All_nodes:
        if len(node) == 4:
            decode += str(Tree.dict_nodes[node]["should"][0])
            Error += (
                abs(Tree.dict_nodes[node]["pred"])
                if sign(Tree.dict_nodes[node]["pred"])
                != Tree.dict_nodes[node]["should"][0]
                else 0
            )
    return Error, decode, Tree


all_trees = []

#Our 8 bits of data before all gambling
input_data = [0, 1, 0, 1, 0, 1, 0, 1]
full_input = [0] * 6 + input_data[:2] + [0, 0] + input_data[2:]
encoded = np.array(BPSK(from_8_to_16(input_data)))

#noise parametres
mu = 0.0 
sigma = 0.2 
noise = np.random.normal(mu, sigma, 16)
noisy_vector = encoded + noise

Example_data = [round(i, 2) for i in noisy_vector]  #Data after padding, encoding, BPSK and noise

#Errors making
Example_data[2] *= -1
Example_data[5] *= -1
Example_data[7] *= -1
#Example_data[10] *= -1

print("BPSK data: ")
print(*encoded, sep="\t")
print("Our data after adding noise: ")
for i in Example_data:
    print(i, end="\t")
print()

answer = []

Tree1 = Tree(Example_data)
Tree1.cur_node = ""
all_trees.append(Tree1)
index = 0
while index < len(all_trees):
    curr_Tree = all_trees[index]
    curr_Tree.cur_node = curr_Tree.nodes.popleft()
    node = curr_Tree.cur_node
    if node == "Finish":
        metrics, code, tree = calculate_metrics(curr_Tree)
        answer.append((metrics, code, tree))
        index += 1

    else:
        if node[-1] == "L":
            node1 = curr_Tree.cur_node = curr_Tree.Tree_level_L()
        else:
            node1 = curr_Tree.cur_node = curr_Tree.Tree_level_R()


sorted_answer = sorted(answer, key=lambda x: x[0])

best = sorted_answer[0][1]


# Определение цветовых кодов ANSI
GREEN = "\033[92m"
RED = "\033[91m"
RESET = "\033[0m"

# Вывод строки "Answer of decoder: " с раскраской элементов
print("Answer of decoder: ", end="\n")
for i in range(len(best)):
    if int(best[i]) == full_input[i]:
        colored_element = GREEN + str(best[i]) + RESET
    else:
        colored_element = RED + str(best[i]) + RESET
    # Добавление табуляции, кроме последнего элемента
    if i < len(best) - 1:
        print(colored_element, end="\t")
    else:
        print(colored_element)

# Вывод строки "Input data: " с раскраской элементов
print("Input data: ", end="\n")
for i in range(len(full_input)):
    if int(best[i]) == full_input[i]:
        colored_element = GREEN + str(full_input[i]) + RESET
    else:
        colored_element = RED + str(full_input[i]) + RESET
    
    # Добавление табуляции, кроме последнего элемента
    if i < len(full_input) - 1:
        print(colored_element, end="\t")
    else:
        print(colored_element)

print("Metrics: ", round(sorted_answer[0][0], 2))