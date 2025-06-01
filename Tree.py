import math
import numpy as np
import random
import copy
import collections

PATH_list = ["L", "LL", "LLL", "LLLL",
        "LLLR", "LLL", "LLR", "LLRL", 
        "LLRR", "LLR", "LL", "LR",
        "LRL", "LRLL", "LRLR", "LRL", 
        "LRR", "LRRL", "LRRR", "LRR",
        "LR", "L", "R", "RL", 
        "RLL", "RLLL", "RLLR", "RLL",
        "RLR", "RLRL", "RLRR", "RLR", 
        "RL", "RR", "RRL", "RRLL",
        "RRLR", "RRL", "RRR", "RRRL", 
        "RRRR", "RRR", "RR", "R",
        "Finish"]
def concatinate(x: list, y: list):
    tmp = [x[i] ^ y[i] for i in range(len(x))]
    for i in y:
        tmp.append(i)
    return tmp

def sign(x):
    return 1 if x >= 0 else -1

def L_func(x, y):
    return (x * y) // abs(x * y) * min(abs(x), abs(y))

def R_func(x, y, b):
    return x + y if b == 0 else y - x

def gaussian(m: list):
    return [random.gauss(0, 2) * i for i in range(len(m))]

def BPSK(m: list):
    a = [1 if x == 0 else 0 for x in m]
    return a

def error(m: list, number=2):
    pass

class Tree:
    dict_nodes = {}
    nodes = collections.deque(PATH_list)
    cur_node = ""
    proceed = set()
    proceed.add("ALL")

    def __init__(self, Elements: list[float]):
        self.dict_nodes = {}
        self.dict_nodes["ALL"] = (Elements, [])

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
                        self.dict_nodes[self.cur_node + "L"][1 if len(self.cur_node) != 3 else "should"],
                        self.dict_nodes[self.cur_node + "R"][1 if len(self.cur_node) != 3 else "should"]
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
                self.dict_nodes[node][1].append(self.dict_nodes[self.cur_node[:-1] + "L"]["should"])
                self.dict_nodes[node][1].append(self.dict_nodes[self.cur_node[:-1] + "R"]["should"])
                #self.cur_node = node
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
                    Help_Result.append(R_func(tmp[i], tmp[i + length // 2], self.dict_nodes[node + "L" if node != "ALL" else "L"][1][i]))
                self.dict_nodes[self.cur_node] = (Help_Result, [])
                return self.cur_node
            else:
                self.dict_nodes[self.cur_node] = (
                    self.dict_nodes[self.cur_node][0],
                    concatinate(
                        self.dict_nodes[self.cur_node + "L"][1 if len(self.cur_node) != 3 else "should"],
                        self.dict_nodes[self.cur_node + "R"][1 if len(self.cur_node) != 3 else "should"]
                    ),
                )
                return self.cur_node
        

def calculate_metrics(tree):
    return 0, 1, 2




all_trees = []

Example_data = [
    0.6,
    0.7,
    0.79,
    0.54,
    0.4,
    0.8,
    -0.3,
    0.39,
    0.44,
    0.36,
    -0.9,
    -0.13,
    -0.81,
    0.62,
    -0.48,
    -0.6,
]
N = 16
K = 8
L = 4



SAFE_BITS = ["LLLL", "LLLR", "LLRL", "LLRR", "LRLL", "LRLR", "RLLL", "RLLR"]

answer = []

Tree1 = Tree(Example_data)
Tree1.cur_node = ""
all_trees.append(Tree1)
index = 0
while index < len(all_trees):
    curr_Tree = all_trees[index]
    
    curr_Tree.cur_node = curr_Tree.nodes.popleft()
    node = curr_Tree.cur_node
    print(index, node)
    if node == "Finish":
        metrics, code, tree = calculate_metrics(curr_Tree)
        answer.append((metrics, code, tree))
        index += 1
        print(index, len(all_trees))
    
    else:
        if node[-1] == "L":
            node1 = curr_Tree.cur_node = curr_Tree.Tree_level_L()
        else:
            node1 = curr_Tree.cur_node = curr_Tree.Tree_level_R()
        

print(*answer, sep="\n")