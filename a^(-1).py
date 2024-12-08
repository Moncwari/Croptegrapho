# Ищем обратный элемент
# К данному числу в кольце

from tabulate import tabulate

n, a = 2467, 877  # Вписываем модуль кольца (n)
# и число (a), обратное к которому мы хотим найти
y2, y1 = 0, 1  # Тут НИЧО НЕ МЕНЯЕМ
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
if d < 0:
    d1 = d + n_old
else:
    d1 = d
print(str(a_old) + "^(-1) mod " + str(n_old) + " = " + str(d) + " = " + str(d1))
