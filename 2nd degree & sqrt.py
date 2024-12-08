# Извлекаем квадратный корень
# По модулю простого числа
# В том числе решаем сравнение
# Второй степени

from tabulate import tabulate
from sympy import legendre_symbol


def reversed(a_1, a_2):  # В функции НИЧО НЕ МЕНЯЕМ
    n, a = a_1, a_2
    y2, y1 = 0, 1
    table = [["q", "r", "y", "n", "a", "y2", "y1"], ["-", "-", "-", n, a, y2, y1]]
    n_old, a_old = n, a

    while a != 0:
        q = n // a
        r = n % a
        y = y2 - q * y1
        n, a, y1, y2 = a, r, y, y1
        listing = [q, r, y, n, a, y2, y1]
        table.append(listing)

    print(tabulate(table, headers="firstrow", tablefmt="fancy_grid"))
    d = y2
    print(str(a_old) + "^(-1) mod " + str(n_old) + " = " + str(d))
    return d


a, p = 4, 5  # Корень из a по модулю p (p - ПРОСТОЕ)
if legendre_symbol(a, p) == -1:
    print("No roots")
    _ = input()

else:
    if legendre_symbol(-1, p) == -1:
        b = -1
        print("b =", b)
    else:
        b = 2
        while legendre_symbol(b, p) != -1:
            b += 1
        print("b =", b)

    s = 0
    p_1 = p - 1
    while p_1 % 2 == 0:
        p_1 = p_1 // 2
        s += 1
    print("s =", s)
    t = (p - 1) // 2**s
    print("t =", t)

    d = reversed(p, a)
    C_0 = (b**t) % p
    print("C_0 =", C_0)
    r = (a ** ((t + 1) // 2)) % p
    print("r =", r)

    for i in range(1, s):
        d_i = ((r**2 * d) ** (2 ** (s - i - 1))) % p
        print("d_" + str(i) + " = " + str(d_i))
        if d_i % p == -1:
            r *= C_0
            print("Updated: r =", r)
        C_0 = (C_0**2) % p
        print("C_0 =", C_0)

    print("Roots: " + str(r) + ", " + str(-r))
