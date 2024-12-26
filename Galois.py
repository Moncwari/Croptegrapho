from sympy.abc import x
import sympy
import itertools
from sympy import Poly, Pow

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
        polynomial = Poly(polynomial, x).set_modulus(self.p).as_expr()
        #print(polynomial)
        if Poly(polynomial, x, modulus = self.p) and self.n == Poly(polynomial).total_degree():
            self.irreducible_polynomial = polynomial
        else: return(sympy.factor(polynomial), "It is an reducible polynomial") 
    
    def Generation_of_irreducible_polynomials(self):
        Irreducible_polynomials = Galois_field(self.p, self.n + 1)
        for i in range(self.p ** (self.n + 1)):
            polynomial = Irreducible_polynomials.Galois_field[i]
            if polynomial == sympy.factor(polynomial) and Poly(polynomial, x).total_degree() == self.n:
                continue
            else:
                Irreducible_polynomials.Galois_field[i] = 0
        Irreducible_polynomials.Galois_field = set(Irreducible_polynomials.Galois_field)
        Irreducible_polynomials.Galois_field.remove(0)
        Irreducible_polynomials.Galois_field = list(Irreducible_polynomials.Galois_field)
        return Irreducible_polynomials.Galois_field
    
    def Addition_or_multiplication_of_polynomials(self, mode: int, polynomials: list):
        result = sympy.sympify(polynomials[0])
        match mode:
            case 0:
                for i in range(1, len(polynomials)):
                    result += sympy.sympify(polynomials[i])
            case 1: # Умножение
                for i in range(1, len(polynomials)):
                    result = result * sympy.sympify(polynomials[i])
                print(result)
                while Poly(result, x).total_degree() >= self.n:
                    result = (Poly(result, x).div(Poly(self.irreducible_polynomial, x)))[1].as_expr()
        result = Poly(result, x).set_modulus(self.p).as_expr()
        return result

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
                while Poly(res, x).total_degree() >= self.n:
                    res = (Poly(res, x).div(Poly(self.irreducible_polynomial, x)))[1]
                    #print(res)
                    res = (res.clear_denoms()[1]).as_expr()
                    #print(res)
                res = Poly(res, x).set_modulus(self.p).as_expr()
                #print(res, "res")
                #print(degree, "degree")
            if degree == (self.p ** self.n) - 1:
                Generators_of_the_field.append(polynom)
        return Generators_of_the_field

    def Decomposition_by_degrees_of_the_selected_generator(self, generative: str):
        generative = sympy.sympify(generative)
        q = (self.p ** self.n) - 1
        generaive_check = Poly(Pow(generative, q), x).div(Poly(self.irreducible_polynomial, x))[1]
        generaive_check = (generaive_check.clear_denoms()[1]).as_expr()
        generaive_check = Poly(generaive_check, x).set_modulus(self.p).as_expr() 
        if generaive_check == 1:
            Decomposition = []
            for i in range(q):
                Element = Poly(Pow(generative, i), x).div(Poly(self.irreducible_polynomial, x))[1]
                Element = (Element.clear_denoms()[1])
                Element_corrected = Poly([coeff % self.p for coeff in Element.all_coeffs()], x, domain='ZZ')
                Decomposition.append(Element_corrected.as_expr())
            return Decomposition
        else: return "Не образующее"
    
Galois = Galois_field(3, 2)
Galois.Obtaining_an_irreducible_polynomial("2*x**2 - 2*x + 1")
print(Galois.Galois_field)
#print(Galois.Addition_or_multiplication_of_polynomials(1, ["3*x**2 + 3*x", "7*x**2 - 5*x + 7", "x**3"]))
print(Galois.Finding_the_generators_of_the_field())
print(Galois.Decomposition_by_degrees_of_the_selected_generator("x"))
#A = Galois.Generation_of_irreducible_polynomials()
#print(A)
#print(Construction_of_the_Galois_field(2, 3))

#B = (Poly("x**4", x).div(Poly("2*x**2 - 2*x + 1", x)))[0].as_expr()
#print(B)
#B = (Poly(B, x).clear_denoms()[1]).as_expr()
#print(B)