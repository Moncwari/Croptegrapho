# Решаем систему сравнений
# С помощью китайской теоремы об остатках

from tabulate import tabulate

equations = []
N = 1
for i in range(3):  # Пишем число сравнений в системе
    listing = input().split(" ")  # Пишем в консоли,
    # чему равно x и по какому модулю ЧЕРЕЗ ПРОБЕЛ,
    # по 1 уравнению на 1 строке
    listing[0] = int(listing[0])
    listing[1] = int(listing[1])
    N *= listing[1]
    equations.append(listing)
print()
print("N =", N)
print()

x2s = []
for i in range(len(equations)):
    N_i = N // equations[i][1]
    print("N_" + str(i + 1) + " =", N_i)
    b = equations[i][1]
    equations[i][1] = N_i
    x2, x1 = 1, 0  # Тут НИЧО НЕ МЕНЯЕМ
    table = [
        ["q", "r", "x", "N_" + str(i + 1), "n_" + str(i + 1), "x2", "x1"],
        ["-", "-", "-", N_i, b, x2, x1],
    ]
    while b != 0:
        q = N_i // b
        r = N_i % b
        x = x2 - q * x1
        N_i, b, x1, x2 = b, r, x, x1
        listing = [q, r, x, N_i, b, x2, x1]
        table.append(listing)

    x2s.append(x2)
    print(tabulate(table, headers="firstrow", tablefmt="fancy_grid"))
    print()

answer = 0
ans_string = "a = ("
for i in range(len(equations)):
    answer += equations[i][0] * x2s[i] * equations[i][1]
    ans_string += (
        str(equations[i][0]) + "*" + str(x2s[i]) + "*" + str(equations[i][1]) + " + "
    )
ans_string = ans_string[:-3] + ") mod " + str(N) + " = " + str(answer % N)
print(ans_string)
