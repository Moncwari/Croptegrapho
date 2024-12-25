from sympy.abc import x
import sympy
import itertools
from sympy import Poly

class Galois_field:
    
        
    def __init__(self, p, n):          # p - модуль, n - степень неприводимого
        self.p = p
        self.n = n
        self.coefficients = [k for k in range(p)]
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

    def Obtaining_an_irreducible_polynomial(self, polynomial: str):
        polynomial = sympy.sympify(polynomial)
        #print(polynomial)
        polynomial = Poly(polynomial, x).set_modulus(self.p).as_expr()
        print(polynomial)
        if polynomial == sympy.factor(polynomial):
            return(polynomial, "It is an irreducible polynomial")
        else: return(sympy.factor(polynomial), "It is an reducible polynomial") 
    

Galois = Galois_field(2, 3)
print(Galois.Obtaining_an_irreducible_polynomial("x**7"))
#print(Galois.Galois_field)
#print(Construction_of_the_Galois_field(2, 3))