from sympy.abc import x
import sympy
import itertools
from sympy import Poly, Pow
import math

def get_padded_coeffs(polynomials, variable):
  """
  Получает список коэффициентов многочленов, дополненных нулями до максимальной степени.

  Args:
    polynomials: Список объектов Poly.
    variable: Символ переменной многочленов.

  Returns:
    Список списков коэффициентов, дополненных нулями.
  """

  if not polynomials:
     return []

  max_degree = max(sympy.degree(poly) for poly in polynomials)
  coeffs_list = []
  for poly in polynomials:
    coeffs = Poly(poly, x).all_coeffs()
    padded_coeffs = coeffs
    
    # Дополняем нулями, если степень меньше максимальной
    if poly == 0:
        num_missing_zeros = 2
    else:
        num_missing_zeros = max_degree - sympy.degree(poly)
    if num_missing_zeros > 0:
      padded_coeffs = [0] * num_missing_zeros + coeffs

    coeffs_list.append(padded_coeffs)
  return coeffs_list

def Text_to_binary(Text: str) -> str:
    return  "".join(f"{ord(i):08b}" for i in Text)

def binary_to_text(binary_string):
  """Переводит двоичную строку в обычную строку."""
  text = ""
  for i in range(0, len(binary_string), 8):
      binary_char = binary_string[i:i+8]
      #print(f"Binary char: {binary_char}") # Выводим отладочную информацию
      if len(binary_char) == 8:
         try:
             text += chr(int(binary_char, 2))
         except ValueError:
             print(f"Error: Could not convert {binary_char} to integer.")
  return text

def Binary_text_to_blocks(Binary_Text: str, n: int, p: int = 2):
    Length_of_binary_text = len(Binary_Text)
    #print(Length_of_binary_text)
    Check_value = 0
    if Length_of_binary_text % n == 0:
        Number_of_blocks = Length_of_binary_text // n
    else: 
        Check_value = 1
        Number_of_blocks = math.ceil(Length_of_binary_text / n)
    Blocks_with_binary_text = []
    for i in range(Number_of_blocks):
        Blocks_with_binary_text.append("")
    match Check_value:
        case 0:
            for i in range(Number_of_blocks):
                for j in range(n):
                    Blocks_with_binary_text[i] += Binary_Text[j + n * i]
        case 1:
            for i in range(Number_of_blocks - 1):
                for j in range(n):
                    Blocks_with_binary_text[i] += Binary_Text[j + n * i]
            for j in range(n - Length_of_binary_text % n - 1):
                #print(i, j, n, len(Binary_Text))
                Blocks_with_binary_text[Number_of_blocks - 1] += Binary_Text[n * (i + 1) + j]
            Blocks_with_binary_text[Number_of_blocks - 1] += "0" * (n - Length_of_binary_text % n)

    Blocks_with_polynoms = [0] * Number_of_blocks
    for i in range(Number_of_blocks):
        coefficients = [int(coef) for coef in Blocks_with_binary_text[i]]
        Blocks_with_polynoms[i] = Poly(coefficients, x).as_expr()
    #print(Blocks_with_binary_text)
    return Blocks_with_polynoms

class Galois_field:    
    def __init__(self, p, n):          # p - модуль, n - степень неприводимого
        self.p = p
        self.n = n
        self.coefficients = [a for a in range(p)]
        field = list(map(list, itertools.product(self.coefficients, repeat=self.n)))
        degree_field = [] 
        for i in range(self.n):
            degree_field += [x ** i]
        result = []
        for i in range(len(field)):
            result += [0]
            for j in range(len(degree_field)):
                result[i] += field[i][j] * degree_field[j]
        self.Galois_field = result
        self.irreducible_polynomial = sympy.sympify("1")

    def Obtaining_an_irreducible_polynomial(self, polynomial: str):
        polynomial = sympy.sympify(polynomial)
        #print(polynomial)
        polynomial = Poly(polynomial, x).set_modulus(self.p)
        Polynomial_corrected = Poly([coeff % self.p for coeff in polynomial.all_coeffs()], x)
        coefficients = Poly(Polynomial_corrected, x).all_coeffs()
        for j in range(1, self.p ** (self.n)): 
                divisor = Poly(self.Galois_field[j], x, domain = 'ZZ')
                if divisor == 0 or Poly(divisor, x).total_degree() >= self.n:
                    continue
                
                remainder = Poly(polynomial,x, domain='ZZ').div(Poly(divisor, x, domain='ZZ'))[1]
                remainder_corrected = Poly([coeff % self.p for coeff in remainder.all_coeffs()], x)
        #print(Poly(Polynomial_corrected, x))
        if Polynomial_corrected in self.Generation_of_irreducible_polynomials() == False:
            return("It is an reducible polynomial 1")
        for i in range(self.p):
            if ((Polynomial_corrected.subs({x: i})) % self.p + self.p) % self.p == 0:
                return("It is an reducible polynomial 2")
        if Poly(polynomial, x, modulus = self.p) and self.n == Poly(polynomial).total_degree() and remainder_corrected != 0:
            self.irreducible_polynomial = polynomial.as_expr()
            return "The polynomial is irreducible, stored in the field."
        else: return("It is an reducible polynomial 3") 
    
    def Generation_of_irreducible_polynomials(self):
        Irreducible_polynomials = Galois_field(self.p, self.n + 1)
        irreducible_polynomials = []
        for i in range(self.p ** (self.n + 1)):
            polynomial = Poly(Irreducible_polynomials.Galois_field[i], x, domain = 'ZZ') 
            if polynomial == 0 or Poly(polynomial, x).total_degree() != self.n: 
                continue
            is_irreducible = True
            for j in range(1, self.p ** (self.n)): 
                divisor = Poly(Irreducible_polynomials.Galois_field[j], x, domain = 'ZZ')
                if divisor == 0 or Poly(divisor, x).total_degree() >= self.n:
                    continue
                
                remainder = Poly(polynomial,x, domain='ZZ').div(Poly(divisor, x, domain='ZZ'))[1]
                remainder_corrected = Poly([coeff % self.p for coeff in remainder.all_coeffs()], x)
                if remainder_corrected == 0:
                  is_irreducible = False
                  break

            if is_irreducible:
                irreducible_polynomials.append(polynomial.as_expr())
        for j in range(len(irreducible_polynomials)):
            for i in range(self.p):
                if ((sympy.sympify(str(irreducible_polynomials[j])).subs({x: i})) % self.p + self.p) % self.p == 0:
                    irreducible_polynomials[j] = -1

        irreducible_polynomials = set(irreducible_polynomials)
        irreducible_polynomials.remove(-1)
        irreducible_polynomials = list(irreducible_polynomials)
        return irreducible_polynomials
    
    def Addition_or_multiplication_of_polynomials(self, mode: int, polynomials: str):   
        polynomials = polynomials.split("/")
        result = sympy.sympify(polynomials[0])
        match mode:
            case 0: # Сложение
                for i in range(1, len(polynomials)):
                    result += sympy.sympify(polynomials[i])
                    while Poly(result, x).total_degree() >= self.n:
                        result = (Poly(result, x).div(Poly(self.irreducible_polynomial, x)))[1]
                        result = (result.clear_denoms()[1]).as_expr()
            case 1: # Умножение
                for i in range(1, len(polynomials)):
                    result = result * sympy.sympify(polynomials[i])
                #print(result)
                while Poly(result, x).total_degree() >= self.n:
                    result = (Poly(result, x).div(Poly(self.irreducible_polynomial, x)))[1]
                    result = (result.clear_denoms()[1]).as_expr()
        result = Poly(result, x).set_modulus(self.p)
        Result_corrected = Poly([coeff % self.p for coeff in result.all_coeffs()], x, domain='ZZ')
        return Result_corrected.as_expr()

    def Finding_the_generators_of_the_field(self):
        Multiplicative_group = set(self.Galois_field)
        Multiplicative_group.remove(0)
        Multiplicative_group = list(Multiplicative_group)
        #print(Multiplicative_group)
        Generators_of_the_field = []
        for polynom in Multiplicative_group:
            degree = 0
            res = 0
            while res != 1:
                degree += 1
                res = Pow(polynom, degree)
                #print(res)
                while Poly(res, x).total_degree() >= self.n:
                    #print(res, "check")
                    #print(Poly(self.irreducible_polynomial, x), "check_2")
                    res = (Poly(res, x).div(Poly(self.irreducible_polynomial, x)))[1]
                    #print(res.as_expr())
                    res = (res.clear_denoms()[1]).as_expr()
                    #print(res)
                res = Poly(res, x).set_modulus(self.p).as_expr()
                #print(res, "res")
                #print(degree, "degree")
                if degree > (self.p ** self.n) - 1:
                    break
            if degree == (self.p ** self.n) - 1:
                Generators_of_the_field.append(polynom)
        return Generators_of_the_field

    def Decomposition_by_degrees_of_the_selected_generator(self, generative: str):
        generative = sympy.sympify(generative)
        q = (self.p ** self.n) - 1
        #generative_check = Poly(Pow(generative, q), x).div(Poly(self.irreducible_polynomial, x))[1]
        #generative_check = (generative_check.clear_denoms()[1]).as_expr()
        #generative_check = Poly(generative_check, x).set_modulus(self.p).as_expr() 
        if Poly(generative, x) in self.Finding_the_generators_of_the_field() == True:
            Decomposition = []
            for i in range(q):
                Element = Poly(Pow(generative, i), x).div(Poly(self.irreducible_polynomial, x))[1]
                Element = (Element.clear_denoms()[1])
                Element_corrected = Poly([coeff % self.p for coeff in Element.all_coeffs()], x, domain='ZZ')
                Decomposition.append(Element_corrected.as_expr())
            return Decomposition
        else: return "Не образующее"

def Affine_encryption(Open_text: str, Galois: Galois_field, a, b) -> str:
    #print(Galois.Galois_field)
    #if Poly(a, x).as_expr() or Poly(b, x).as_expr() in Galois.Galois_field == False:
    #    return("a или b не принадлежат полю Галуа.")
    Binary_open_text = Text_to_binary(Open_text)
    #print(Binary_open_text)
    Text_for_encryption = Binary_text_to_blocks(Binary_open_text, Galois.n) # По сути х(х1, х2, ...)
    #print(Text_for_encryption)
    #print(Text_for_encryption)
    Encryption_text = ""
    Blocks_with_encryption_text = []
    for i in range(len(Text_for_encryption)):
        Element = a * Text_for_encryption[i] + b
        #print(Element)
        Element = Poly(Element, x).div(Poly(Galois.irreducible_polynomial, x))[1]
        Element = (Element.clear_denoms()[1])
        Element_corrected = Poly([coeff % Galois.p for coeff in Element.all_coeffs()], x)
        Blocks_with_encryption_text.append(Element_corrected.as_expr())
    for i in range(len(Blocks_with_encryption_text)):
        Encryption_text += str(Blocks_with_encryption_text[i]) + "/"
    return Encryption_text

def Affine_decryption(Cipher_text: str, Galois: Galois_field, a, b):
    Open_text = ""
    #print(Cipher_text)
    Blocks_with_cipher_text = Cipher_text.split("/") # x, правда пока не преобразованы в полиномы
    #print(Blocks_with_cipher_text)
    for a_1 in Galois.Galois_field:
        Element = a_1 * a
        #print(Element, a_1)
        Element = Poly(Element, x).div(Poly(Galois.irreducible_polynomial, x))[1]
        Element = (Element.clear_denoms()[1])
        Element_corrected = Poly([coeff % Galois.p for coeff in Element.all_coeffs()], x).as_expr()
        if Element_corrected == 1:
            break
    #print(a_1)
    Blocks_with_decryption_text = []
    for i in range(len(Blocks_with_cipher_text)-1):
        #print(Blocks_with_cipher_text[i])
        Element = (sympy.sympify(Blocks_with_cipher_text[i]) - b) * a_1
        #print(Element)
        Element = Poly(Element, x).div(Poly(Galois.irreducible_polynomial, x))[1]
        Element = (Element.clear_denoms()[1])
        Element_corrected = Poly([coeff % Galois.p for coeff in Element.all_coeffs()], x)
        Blocks_with_decryption_text.append(Element_corrected.as_expr())
    Blocks_with_coefficients = get_padded_coeffs(Blocks_with_decryption_text, x)
    #print(Blocks_with_coefficients)
    #for i in range(len(Blocks_with_decryption_text)):
    #    Blocks_with_coefficients.append(Poly(Blocks_with_decryption_text[i], x).all_coeffs())
    #    print(Blocks_with_decryption_text[i], Poly(Blocks_with_decryption_text[i], x).all_coeffs())
    Blocks_with_binary_text = []
    for i in range(len(Blocks_with_coefficients)):
        Element = ""
        for j in range(Galois.n):
            Element += str(Blocks_with_coefficients[i][j])
        Blocks_with_binary_text.append(Element)
    #print(Blocks_with_binary_text)
    for i in range(len(Blocks_with_binary_text)):
        Open_text += Blocks_with_binary_text[i]
    #print(Open_text)
    while len(Open_text) % 8 != 0:
        Open_text = Open_text[0: len(Open_text) - 1]
    #print(Open_text)
    return binary_to_text(Open_text)

#Galois = Galois_field(2, 3)
#Galois.Obtaining_an_irreducible_polynomial("x**3 + x**2 + x + 1")
#print(Galois.Obtaining_an_irreducible_polynomial("x**3"))
#print(Galois.Galois_field)
#print(Affine_encryption("hello world", Galois, x, 1))
#print(Affine_decryption(Affine_encryption("hello world", Galois, x, 1), Galois, x, 1))
#print(Galois.Addition_or_multiplication_of_polynomials(1, ["3*x**2 + 3*x", "7*x**2 - 5*x + 7", "x**3"]))
#print(Galois.Finding_the_generators_of_the_field())
#print(Galois.Decomposition_by_degrees_of_the_selected_generator("x"))
#A = Galois.Generation_of_irreducible_polynomials()
#print(A)
#print(Construction_of_the_Galois_field(2, 3))

#B = (Poly("x**4", x).div(Poly("2*x**2 - 2*x + 1", x)))[0].as_expr()
#print(B)
#B = (Poly(B, x).clear_denoms()[1]).as_expr()
#print(B)

#Result = Binary_text_to_blocks(Text_to_binary("hello world"), 3)
#for i in range(len(Result)):
#    print(Result[i])
Galois = Galois_field(3, 2)
Galois.Obtaining_an_irreducible_polynomial("2*x**2 - 2*x + 1")
#print("Образующие заданного поля:")
#print(Galois.Finding_the_generators_of_the_field())
#print(Galois.Obtaining_an_irreducible_polynomial("2*x**2 - 2*x + 1"))
#print("Введите неприводимый многочлен, умножение - '*', возведение в степень - '**'. Сложение и вычитание стандартные, записываются через пробел. Переменная - х:")
#print("Введите 1, если хотите ввести неприводимый многочлен самостоятельно. Введите 2, если хотите получить список подходящих неприводимых многочленов:")
#print("Введите 0, если хотите сложить многочлены, введите 1, если хотите их умножить:")
#print("Введите обраующее, по степеням которого хотите выполнить разложение:")
print("Введите значения p, n и неприводимый многочлен. Всё в одну строчку через '/':")
input()
print("Введите значения a и b, они должны принадлежать полю Галуа:")
input()
print("Введите шифртекст, который хотите расшифровать:")
input()
#print(Affine_encryption("hello world", Galois, 1, 2))
print(Affine_decryption(Affine_encryption("hello world", Galois, 1, 2), Galois, 1, 2))
#print(Galois.Decomposition_by_degrees_of_the_selected_generator("x + 1"))
#print("Введите многочлены, которые хотите умножить. умножение - '*', возведение в степень - '**'. Сложение и вычитание стандартные, записываются через пробел. Переменная - х. Многочлены разделять '/'.")
#input()
#print(Galois.Addition_or_multiplication_of_polynomials(1, "2*x + 1/x"))
#print(Galois.Generation_of_irreducible_polynomials())
#print(Galois.Obtaining_an_irreducible_polynomial("2*x**2 + 2*x + 2"))
#print("Введите значения p и n в одну строчку через пробел:")
#k = input()
#print(Galois.Galois_field)