from GaloisMoncwari import *

#1
'''
p, n = map(int, input("Enter the modulus and degree:\n").split())
Galois = Galois_field(p, n)
print(Galois.Galois_field)
'''
#2
'''
p, n = 2, 8
Galois = Galois_field(p, n)
poly = input("Enter the irreducible polynomial:\n")
Galois.Obtaining_an_irreducible_polynomial(poly)
'''
#3
'''
p, n = map(int, input("Enter the modulus and degree:\n").split())
Galois = Galois_field(p, n)
A = Galois.Generation_of_irreducible_polynomials()
print("The irreducible polynomials are:\n", A, sep='')
'''
#4
'''
p, n = map(int, input("Enter the modulus and degree:\n").split())
Galois = Galois_field(p, n)
#print(Galois.Generation_of_irreducible_polynomials())
Galois.Obtaining_an_irreducible_polynomial(input("Enter the irreducible polynomial:\n"))
m = int(input("Enter 0 for addition and 1 for multiplication:\n"))
st = input("Enter the polynomials. Polynomials must be separated by /:\n")
print("The result is:\n", Galois.Addition_or_multiplication_of_polynomials(m, list(st.split('/'))), sep='')
'''
#5
'''
p, n = map(int, input("Enter the modulus and degree:\n").split())
Galois = Galois_field(p, n)
Galois.Obtaining_an_irreducible_polynomial("x**3 + x + 1")
print("Generators of the field are:\n", Galois.Finding_the_generators_of_the_field(), sep = '')
'''
#6
'''
p, n = map(int, input("Enter the modulus and degree:\n").split())
Galois = Galois_field(p, n)
Galois.Obtaining_an_irreducible_polynomial("x**2 + 1")
print("The decomposition is:\n", Galois.Decomposition_by_degrees_of_the_selected_generator(input("Enter the generator:\n")), sep = '')
'''