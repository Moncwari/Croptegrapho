import itertools
from sympy import symbols, Poly, factorint, simplify, factor, div
import sympy
from sympy.abc import x

p = 2
n = 3
elements = [k for k in range(p)]  # список элементов поля GF(p)
field = list(map(list, itertools.product(elements, repeat=n)))

for i in range(len(field)):
    print(field[i])

degree_field = [] * n
for i in range(n):
    degree_field += [x ** i]

#degree_field = list(map(list, itertools.product(degree_field, repeat=n)))

#print(degree_field)

#field = list(map(list, itertools.product(elements, degree_field)))
result = []
for i in range(len(field)):
    result += [0]
    for j in range(len(degree_field)):
        result[i] += field[i][j] * degree_field[j]
if result[1] == factor(result[1]):
    print("yeahhh")
print(result)
#for i in range(len(field)):
#    print(field[i])


