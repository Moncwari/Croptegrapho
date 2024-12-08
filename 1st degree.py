# Решаем сравнение первой степени

from tabulate import tabulate

a, b, c = 3, 2, 5  # Пишем по порядку:
# кэф при иксе (a), чему равно (b),
# по какому модулю (c)
x2, x1, y2, y1 = 1, 0, 0, 1  # Тут НИЧО НЕ МЕНЯЕМ
table = [
    ["q", "r", "x", "y", "a", "b", "x2", "x1", "y2", "y1"],
    ["-", "-", "-", "-", c, a, x2, x1, y2, y1],
]
c_old, a_old = c, a

while a != 0:
    q = c // a
    r = c % a
    x = x2 - q * x1
    y = y2 - q * y1
    c, a, x1, x2, y1, y2 = a, r, x, x1, y, y1
    listing = [q, r, x, y, c, a, x2, x1, y2, y1]
    table.append(listing)

d, x, y = c, x2, y2
print(tabulate(table, headers="firstrow", tablefmt="fancy_grid"))
print()

if d == 1:
    print("Nashe sravnenie nelzya sokratit")
    print(str(a_old) + "x mod " + str(c_old) + " = " + str(b))
    print(str(a_old) + "^(-1) mod " + str(c_old) + " = " + str(y))
    print("x mod " + str(c_old) + " = " + str((y * b) % c_old))
else:
    print("Sokratim nashe sravnenie")
    new_a, new_b, new_c = a_old // d, b // d, c_old // d
    print(str(new_a) + "x mod " + str(new_c) + " = " + str(new_b))
    for i in range(1, new_c):
        if (new_a * i) % new_c == 1:
            print(str(new_a) + "^(-1) mod " + str(new_c) + " = " + str(i))
            x = (i * new_b) % new_c
            break
    k = 0
    while x + k * new_c < c_old:
        print("x mod " + str(c_old) + " = " + str(x + k * new_c))
        k += 1
