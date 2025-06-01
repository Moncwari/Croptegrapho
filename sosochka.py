import math


def signed_min(x, y):
    product_sign = 1 if x * y >= 0 else -1  # Определяем знак произведения
    minimum = min(abs(x), abs(y))  # Находим минимум модулей
    return product_sign * minimum


def L(y):
    k = len(y) // 2
    res = []
    for i in range(k):
        res.append(round(signed_min(y[i], y[i + k]), 4))
    print(res)
    return res


def R(xy, b):
    k = len(xy) // 2
    res = []
    for i in range(k):
        if b[i] == 0:
            res.append(round(xy[i + k] + xy[i], 4))
        elif b[i] == 1:
            res.append(round(xy[i + k] - xy[i], 4))
    print(res)
    return res


def add_lists_elementwise(list1, list2):
    if len(list1) != len(list2):
        print("Ошибка: Списки должны быть одинаковой длины для поэлементного сложения.")
        return None

    result_list = []
    for i in range(len(list1)):
        result_list.append((list1[i] + list2[i]) % 2)

    return result_list


u = "0000001000110010"
# u = '0000000100010101'
print(len(u))
l = [[] for _ in range(int(math.sqrt(len(u))))]
print(l)
k = int(math.sqrt(len(u)))
print(k)
for _ in range(int(math.sqrt(len(u)))):
    k -= 1
    if k < 0:
        break
    elif k == int(math.sqrt(len(u))) - 1:
        for i in range(0, len(u), 2):
            l[k].append([(int(u[i]) + int(u[i + 1])) % 2, int(u[i + 1])])
    else:
        for j in range(0, len(l[k + 1]), 2):
            l[k].append(
                add_lists_elementwise(l[k + 1][j], l[k + 1][j + 1]) + l[k + 1][j + 1]
            )
    print(l)
print(f"result:{l[0][0]}")
"""
F = [1, 2, 3, 4, 5, 6, 9, 10]
y = [-1.7, -0.63, -1.3, 1.1, 0.39, -0.91, -1.22, -0.77]
print(len(y), y)
ys = [[] for _ in range(int(math.sqrt(len(y))))]
y1 = L(y)
y2 = L(y1)
y3 = L(y2)
yr1 = R(y2, [0])
yr2 = R(y1, [0, 0])
y4 = L(yr2)
yr3 = R(yr2, [0])"""

# y = [0.9, -1.5, -0.37, -1.1, 1.34, 1.22, 1.75, -1.32,
# -0.84, -0.92, -1.35, -0.52, -0.78, 1.62, -0.23, 1.01]
# Anna 34
y = [
    0.85,
    -0.9,
    -2.3,
    0.21,
    -0.2,
    -0.71,
    0.92,
    -0.34,
    -1.3,
    -0.98,
    1.75,
    1.23,
    1.07,
    0.81,
    -1.3,
    2,
]
# Anna 4
print("----------------------------------")
y1 = L(y)
y2 = L(y1)
y3 = L(y2)
y4 = L(y3)
yr1 = R(y3, [0])
yr2 = R(y2, [0, 0])
y5 = L(yr2)
yr3 = R(yr2, [0])
yr4 = R(y1, [0, 0, 0, 0])
y6 = L(yr4)
y7 = L(y6)
yr5 = R(y6, [0, 0])
yr6 = R(yr4, [0, 0])
y8 = L(yr6)
yr7_1, yr7_2 = R(yr6, [0]), R(yr6, [1])
# Правая часть
yr8_1, yr8_2, yr8_3, yr8_4 = (
    R(y, [0] * 8),
    R(y, [1] * 8),
    R(y, [1, 0] * 4),
    R(y, [0, 1] * 4),
)
y9_1, y9_2, y9_3, y9_4 = L(yr8_1), L(yr8_2), L(yr8_3), L(yr8_4)
y10_1, y10_2, y10_3, y10_4 = (
    L(
        y9_1,
    ),
    L(y9_2),
    L(y9_3),
    L(y9_4),
)
print("---")
y11_1, y11_2, y11_3, y11_4 = L(y10_1), L(y10_2), L(y10_3), L(y10_4)
yr9_1, yr9_2, yr9_3, yr9_4 = R(y10_1, [0]), R(y10_2, [0]), R(y10_3, [0]), R(y10_4, [0])
yr10_1, yr10_2, yr10_3, yr10_4 = (
    R(y9_1, [0] * 2),
    R(y9_2, [0] * 2),
    R(y9_3, [0] * 2),
    R(y9_4, [0] * 2),
)
y12_1, y12_2, y12_3, y12_4 = L(yr10_1), L(yr10_2), L(yr10_3), L(yr10_4)
# Важно!!!!
y13_1, y13_2, y13_3, y13_4 = (
    R(yr10_2, [0]),
    R(yr10_2, [1]),
    R(yr10_3, [0]),
    R(yr10_3, [1]),
)
print("---")
yr11_1, yr11_2, yr11_3, yr11_4 = (
    R(yr8_1, [1] * 4),
    R(yr8_2, [0] * 4),
    R(yr8_3, [0] * 4),
    R(yr8_3, [0, 1] * 2),
)
y14_1, y14_2, y14_3, y14_4 = L(yr11_1), L(yr11_2), L(yr11_3), L(yr11_4)
y15_1, y15_2, y15_3, y15_4 = L(y14_1), L(y14_2), L(y14_3), L(y14_4)
yr12_1, yr12_2, yr12_3, yr12_4 = (
    R(y14_1, [0]),
    R(y14_2, [0]),
    R(y14_4, [0]),
    R(y14_4, [1]),
)
yr13_1, yr13_2, yr13_3, yr13_4 = (
    R(yr11_1, [0] * 2),
    R(yr11_2, [0] * 2),
    R(yr11_4, [0] * 2),
    R(yr11_4, [0, 1]),
)
y16_1, y16_2, y16_3, y16_4 = L(yr13_1), L(yr13_2), L(yr13_3), L(yr13_4)
yr14_1, yr14_2, yr14_3, yr14_4 = (
    R(yr13_1, [0]),
    R(yr13_2, [1]),
    R(yr13_4, [1]),
    R(yr13_4, [0]),
)
