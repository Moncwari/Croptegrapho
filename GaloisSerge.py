from sympy.abc import x
import sympy
from itertools import product
from sympy import Poly, Pow, sympify, degree
import math

class GaloisField:
    def __init__(self, p, n):  # (p-1) - макс.коэф, n - макс.степ
        self.p = p
        self.n = n

    def field_creator(self):
        coeffs = product(range(self.p), repeat=self.n)
        polys = [Poly.from_list(coeff, gens=x) for coeff in coeffs]
        return polys

    def is_irreducible(self, poly):
        poly = Poly(sympify(poly), x).set_modulus(self.p)
        if degree(poly) != self.n and poly.LC() >= p:
            return False
        else:
            fact = poly.factor_list()[1]
            if len(fact) == 1:
                if fact[0][1] == 1:
                    return True
                else:
                    return False
            else:
                return False

    def gen_irreducible(self):
        coeffs = list(product(range(self.p), repeat=(self.n+1)))
        polys = []
        for i in range (len(coeffs)):
            if coeffs[i][0] != 0:
                poly = Poly.from_list(coeffs[i], gens=x)
                if self.is_irreducible(poly):
                    polys.append(poly)
        return polys

    def math_oper(self, first, second, oper, poly1):
        self.first = Poly(sympify(first), x).set_modulus(self.p)
        self.second = Poly(sympify(second), x).set_modulus(self.p)
        self.oper = oper
        self.poly1 = Poly(sympify(poly1), x).set_modulus(self.p)
        if self.oper == "+":
            ans = Poly(self.first + self.second, x).set_modulus(self.p)
            ans = Poly([coeff % self.p for coeff in ans.all_coeffs()], x)
            return ans
        if self.oper == "-":
            ans = Poly(self.first - self.second, x).set_modulus(self.p)
            ans = Poly([coeff % self.p for coeff in ans.all_coeffs()], x)
            return ans
        if self.oper == "*":
            ans = Poly(self.first * self.second, x)
            ans = Poly(ans, x).div(Poly(self.poly1, x))[1]
            ans = Poly([coeff % self.p for coeff in ans.all_coeffs()], x)
            return ans

    def mult_gens(self, poly2):
        self.poly2 = Poly(sympify(poly2), x).set_modulus(self.p)
        elems = self.field_creator()[1:]
        gens = []
        for elem in elems:
            res = 0
            degree = 0
            while res != 1:
                degree += 1
                res = Poly(elem.pow(degree), x).set_modulus(self.p)
                res = Poly(res, x).div(Poly(self.poly2, x))[1]
                if res == 1 and degree == p**n-1:
                    gens.append(elem)
        return gens

    def power_by_gen(self, poly3, gen):
        self.poly3 = Poly(sympify(poly3), x).set_modulus(self.p)
        self.gen = Poly(sympify(gen), x).set_modulus(self.p)
        powers = []
        for i in range (1, p**n):
            power = Poly(self.gen.pow(i), x).set_modulus(self.p)
            power = Poly(power, x).div(Poly(self.poly3, x))[1]
            power = Poly([coeff % self.p for coeff in power.all_coeffs()], x)
            powers.append(power)
        return powers

def encrypt_decrypt(string, alpha, beta, poly, key):
    field = GaloisField(3, 3)
    elems = field.field_creator()
    alphabet = "abcdefghijklmnopqrstuvwxyz_"
    to_Galois = []
    for letter in string:
        index = alphabet.find(letter)
        to_Galois.append(elems[index])

    if key == 1:
        for i in range (len(to_Galois)):
            to_Galois[i] = field.math_oper(alpha, to_Galois[i], "*", poly)
            to_Galois[i] = field.math_oper(to_Galois[i], beta, "+", poly)
    if key == 2:
        rev_alpha = 0
        mult = elems[1:]
        for elem in mult:
            if field.math_oper(alpha, elem, "*", poly) == 1:
                rev_alpha = elem
        for i in range (len(to_Galois)):
            to_Galois[i] = field.math_oper(to_Galois[i], beta, "-", poly)
            to_Galois[i] = field.math_oper(to_Galois[i], rev_alpha, "*", poly)            

    from_Galois = ""
    for elem in to_Galois:
        index = elems.index(elem)
        from_Galois += alphabet[index]
    return from_Galois

print("""Введите 1, если хотите поработать с полем Галуа
      или 2, если хотите произвести шифрование""")
a = int(input())
if a == 1:
    print("""Введите натуральные значения p и n через пробел,
          p - строго простое""")
    string = input().split(" ")
    p, n = int(string[0]), int(string[1])
    field = GaloisField(p, n)
    print("""Введите 1, если хотите построить поле Галуа
          по значениям p и n""")
    print("""Введите 2, если хотите поработать с неприводимыми
          многочленами в поле Галуа""")
    print("""Введите 3, если хотите сложить или умножить
          многочлены в поле Галуа""")
    print("""Введите 4, если хотите поработать с образующими
          элементами мультипликативной группы поля Галуа""")
    b = int(input())
    if b == 1:
        print("Элементы поля Галуа F(" + str(p) + "**" + str(n) + ")")
        elems = field.field_creator()
        for elem in elems:
            print(elem.as_expr(), sep="")
    if b == 2:
        print("""Введите 1, если хотите проверить ваш многочлен
              на неприводимость""")
        print("""Введите 2, если хотите получить список всех неприводимых
              многочленов""")
        c = int(input())
        if c == 1:
            print("""Введите неприводимый многочлен степени n
                  с коэффицентами не более p-1""")
            string = input()
            if field.is_irreducible(string):
                print("Многочлен действительно неприводим")
            else:
                print("Что-то пошло не так...")
        if c == 2:
            irrs = field.gen_irreducible()
            print("Все неприводимые многочлены степени n")
            for irr in irrs:
                print(irr.as_expr(), sep="")
    if b == 3:
        print("""Введите 2 элемента поля Галуа, по 1 в строке,
              в строке между ними - знак арифметической операции,
              затем - неприводимый многочлен степени n с коэффицентами
              не более p-1""")
        first = input()
        oper = input()
        second = input()
        poly = input()
        ans = field.math_oper(first, second, oper, poly)
        print(ans.as_expr())
    if b == 4:
        print("""Введите неприводимый многочлен степени n
              с коэффицентами не более p-1""")
        poly = input()
        elems = field.mult_gens(poly)
        print("Образующие элементы мультипликативной группы")
        for elem in elems:
            print(elem.as_expr(), sep="")
        print("Выберите один из образующих элементов")
        gen = input()
        print("Разложение по степеням выбранного образующего")
        powers = field.power_by_gen(poly, gen)
        for i in range(1, p**n):
            print("(" + str(gen) + ")" + "**" + str(i) + " = "
                  + str(powers[i-1].as_expr()))
if a == 2:
    print("""Шифрование производится в алфавите
          abcdefghijklmnopqrstuvwxyz_""")
    print("Используется поле Галуа F(3**3)")
    field = GaloisField(3, 3)
    elems = field.field_creator()
    for elem in elems:
        print(elem.as_expr(), sep="")
    print("""Введите alpha и beta, принадлежащие F(3**3),
          alpha != 0, по 1 в строке""")
    alpha = input()
    beta = input()
    print("Список неприводимых многочленов степени 3")
    irrs = field.gen_irreducible()
    for irr in irrs:
        print(irr.as_expr(), sep="")
    print("""Введите неприводимый многочлен степени 3
              с коэффицентами не более 2""")
    poly = input()
    print("Введите ваш текст")
    text = input()
    print("""Введите 1, если хотите зашифровать его,
          или 2, если расшифровать""")
    b = int(input())
    print(encrypt_decrypt(text, alpha, beta, poly, b))
