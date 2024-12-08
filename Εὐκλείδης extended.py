# Ищем наибольший общий делитель 2 чисел
# И их линейную комбинацию

from tabulate import tabulate

a, b = 3931, 1148  # У чисел a и b необходимо найти НОД
x2, x1, y2, y1 = 1, 0, 0, 1  # Тут НИЧО НЕ МЕНЯЕМ
table = [
    ["q", "r", "x", "y", "a", "b", "x2", "x1", "y2", "y1"],
    ["-", "-", "-", "-", a, b, x2, x1, y2, y1],
]
a_old, b_old = a, b

while b != 0:
    q = a // b
    r = a % b
    x = x2 - q * x1
    y = y2 - q * y1
    a, b, x1, x2, y1, y2 = b, r, x, x1, y, y1
    listing = [q, r, x, y, a, b, x2, x1, y2, y1]
    table.append(listing)

print(tabulate(table, headers="firstrow", tablefmt="fancy_grid"))
d, x, y = a, x2, y2
print("NOD raven", d)
print(str(x) + "*" + str(a_old) + " + " + str(y) + "*" + str(b_old) + " = " + str(d))
