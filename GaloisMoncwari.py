from sympy.abc import x
import sympy
from itertools import product
from sympy import Poly, Pow
from math import ceil
import sys


def get_padded_coeffs(polynomials, n):
    """
    Gets a list of coefficients of polynomials, padded with zeros to the maximum degree.

    Args:
        polynomials: A list of Poly objects.
        variable: The variable symbol of the polynomials.

    Returns:
        A list of lists of coefficients, padded with zeros.
    """
    if not polynomials:
        return []

    coeffs_list = []
    for poly in polynomials:
        coeffs = Poly(poly, x).all_coeffs()
        padded_coeffs = coeffs

        # Pad with zeros if the degree is less than the maximum
        if poly == 0:
            num_missing_zeros = n - 1
        else:
            num_missing_zeros = n - sympy.degree(poly) - 1
        if num_missing_zeros > 0:
            padded_coeffs = [0] * num_missing_zeros + coeffs

        coeffs_list.append(padded_coeffs)
    return coeffs_list


def Text_to_binary(Text: str) -> str:
    """
    Converts a given text string to a binary string.

    Args:
        Text: The input text to convert.

    Returns:
        A string representing the binary equivalent of the input text.
    """
    # Convert each character in the text to its ASCII binary representation
    return "".join(format(ord(i), "08b") for i in Text)


def Binary_to_text(binary_string):
    """
    Translates a binary string to a regular string.

    :param binary_string: A binary string.
    :return: A regular string.
    """
    try:
        # Convert each 8-bit binary number to a character
        return "".join(
            chr(int(binary_string[i : i + 8], 2))
            for i in range(0, len(binary_string), 8)
        )
    except:
        # Return an error message if there's a problem
        return f"You can't convert {binary_string} to text"


def Binary_text_to_blocks(Binary_Text: str, n: int) -> list:
    """
    Breaks down a binary string into blocks of length n.

    Args:
        Binary_Text (str): A binary string.
        n (int): The length of each block.

    Returns:
        list: A list of blocks, each of which is a string of length n.

    """
    Length_of_binary_text = len(Binary_Text)
    Check_value = 0
    if Length_of_binary_text % n == 0:
        Number_of_blocks = Length_of_binary_text // n
    else:
        Check_value = 1
        Number_of_blocks = ceil(Length_of_binary_text / n)
    Blocks_with_binary_text = [""] * Number_of_blocks
    match Check_value:
        case 0:
            # If the length of the binary text is divisible by n, then simply
            # split the string into blocks of length n
            for i in range(Number_of_blocks):
                for j in range(n):
                    Blocks_with_binary_text[i] += Binary_Text[j + n * i]
        case 1:
            # If the length of the binary text is not divisible by n, then
            # split the string into blocks of length n and pad the last block
            # with zeros to make it a length of n
            for i in range(Number_of_blocks - 1):
                for j in range(n):
                    Blocks_with_binary_text[i] += Binary_Text[j + n * i]
            for j in range(n - Length_of_binary_text % n - 1):
                Blocks_with_binary_text[Number_of_blocks - 1] += Binary_Text[
                    n * (i + 1) + j
                ]
            Blocks_with_binary_text[Number_of_blocks - 1] += "0" * (
                n - Length_of_binary_text % n
            )

    Blocks_with_polynoms = [0] * Number_of_blocks
    for i in range(Number_of_blocks):
        coefficients = [int(coef) for coef in Blocks_with_binary_text[i]]
        Blocks_with_polynoms[i] = Poly(coefficients, x).as_expr()
    return Blocks_with_polynoms


class Galois_field:
    def __init__(self, p, n):
        """
        Initializes a Galois field of degree n over the modulus p.

        Args:
            p: The modulus of the Galois field
            n: The degree of the Galois field

        Attributes:
            p: The modulus of the Galois field
            n: The degree of the Galois field
            coefficients: A list of the coefficients of the elements of the Galois field
            Galois_field: A list of the elements of the Galois field
            irreducible_polynomial: The irreducible polynomial of the Galois field
        """
        self.p = p
        self.n = n
        self.coefficients = [a for a in range(p)]
        self.Galois_field = [
            sum(coeff * x**i for i, coeff in enumerate(field))
            for field in product(self.coefficients, repeat=self.n)
        ]
        self.irreducible_polynomial = sympy.sympify("1")

    def Obtaining_an_irreducible_polynomial(self, polynomial: str):
        polynomial = sympy.sympify(polynomial)
        fl = True
        polynomial = Poly(polynomial, x).set_modulus(self.p)
        Polynomial_corrected = Poly(
            [coeff % self.p for coeff in polynomial.all_coeffs()], x
        )
        for j in range(1, self.p ** (self.n)):
            divisor = Poly(self.Galois_field[j], x, domain="ZZ")
            if divisor == 0 or Poly(divisor, x).total_degree() >= self.n:
                continue

            remainder = Poly(polynomial, x, domain="ZZ").div(
                Poly(divisor, x, domain="ZZ")
            )[1]
            remainder_corrected = Poly(
                [coeff % self.p for coeff in remainder.all_coeffs()], x
            )
        if (
            Polynomial_corrected
            in self.Generation_of_irreducible_polynomials()
            == False
        ):
            fl = False
        for i in range(self.p):
            if ((Polynomial_corrected.subs({x: i})) % self.p + self.p) % self.p == 0:
                fl = False
        if (
            Poly(polynomial, x, modulus=self.p)
            and self.n == Poly(polynomial).total_degree()
            and remainder_corrected != 0
        ):
            self.irreducible_polynomial = polynomial.as_expr()
        else:
            fl = False
        if fl:
            print("Polynomial is irreducible")
        else:
            print("Polynomial is reducible")
            sys.exit()

    def Generation_of_irreducible_polynomials(self):
        """
        Finds all the irreducible polynomials of a given degree over the Galois field.

        Returns:
            A list of irreducible polynomials as strings
        """
        # Create a Galois field of the same modulus, but with degree n+1
        Irreducible_polynomials = Galois_field(self.p, self.n + 1)
        irreducible_polynomials = []
        # Iterate over all the elements of the Galois field
        for i in range(self.p ** (self.n + 1)):
            polynomial = Poly(Irreducible_polynomials.Galois_field[i], x, domain="ZZ")
            if polynomial == 0 or Poly(polynomial, x).total_degree() != self.n:
                continue  # Check if the polynomial has the right degree
            is_irreducible = True  # Assume the polynomial is irreducible
            for j in range(1, self.p ** (self.n)):
                divisor = Poly(Irreducible_polynomials.Galois_field[j], x, domain="ZZ")
                if divisor == 0 or Poly(divisor, x).total_degree() >= self.n:
                    continue  # Check if the divisor has the right degree

                remainder = Poly(polynomial, x, domain="ZZ").div(
                    Poly(divisor, x, domain="ZZ")
                )[1]
                remainder_corrected = Poly(
                    [coeff % self.p for coeff in remainder.all_coeffs()], x
                )
                if remainder_corrected == 0:
                    is_irreducible = False
                    break  # If the remainder is zero, the polynomial is not irreducible

            if is_irreducible:
                irreducible_polynomials.append(polynomial.as_expr())

        return irreducible_polynomials

    def Addition_or_multiplication_of_polynomials(self, mode: int, polynomials: list):
        """
        Perform addition or multiplication of polynomials over the Galois field.

        Args:
            mode: 0 for addition and 1 for multiplication
            polynomials: A list of polynomials as strings

        Returns:
            The result of the operation as a polynomial string
        """
        result = sympy.sympify(polynomials[0])
        match mode:
            case 0:  # Addition
                for i in range(1, len(polynomials)):
                    result += sympy.sympify(polynomials[i])
            case 1:  # Multiplication
                for i in range(1, len(polynomials)):
                    result = result * sympy.sympify(polynomials[i])
                # Reduce the result modulo the Galois field
                while Poly(result, x).total_degree() >= self.n:
                    result = (
                        Poly(result, x).div(Poly(self.irreducible_polynomial, x))
                    )[1].as_expr()
        # Set the modulus of the result to the Galois field
        result = Poly(result, x).set_modulus(self.p).as_expr()
        return result

    def Finding_the_generators_of_the_field(self):
        """
        Finds all the generators of the Galois field.

        Returns:
            A list of all the generators of the Galois field.
        """
        Multiplicative_group = set(self.Galois_field)
        Multiplicative_group.remove(0)
        Multiplicative_group = list(Multiplicative_group)
        Generators_of_the_field = []
        for polynom in Multiplicative_group:
            degree = 0
            res = 0
            # Check if the element is a generator of the Galois field
            while res != 1:
                degree += 1
                res = Pow(polynom, degree)
                # Reduce the element modulo the irreducible polynomial
                while Poly(res, x).total_degree() >= self.n:
                    res = (Poly(res, x).div(Poly(self.irreducible_polynomial, x)))[1]
                    res = (res.clear_denoms()[1]).as_expr()
                # Reduce the element modulo the prime
                res = Poly(res, x).set_modulus(self.p).as_expr()
                if degree > (self.p**self.n) - 1:
                    break
            if degree == (self.p**self.n) - 1:
                Generators_of_the_field.append(polynom)
        return Generators_of_the_field

    def Decomposition_by_degrees_of_the_selected_generator(self, generative: str):
        """
        Finds the decomposition of a given generative element of the Galois field.

        Args:
            generative: The generative element as a string

        Returns:
            A list of the decomposition of the generative element by degrees
        """
        generative = sympy.sympify(generative)
        q = (self.p**self.n) - 1

        # Check if the element is a generator of the Galois field
        generaive_check = Poly(Pow(generative, q), x).div(
            Poly(self.irreducible_polynomial, x)
        )[1]
        generaive_check = (generaive_check.clear_denoms()[1]).as_expr()
        generaive_check = Poly(generaive_check, x).set_modulus(self.p).as_expr()
        if generaive_check == 1:
            Decomposition = []
            for i in range(q):
                # Find the decomposition of the element by degrees
                Element = Poly(Pow(generative, i), x).div(
                    Poly(self.irreducible_polynomial, x)
                )[1]
                Element = Element.clear_denoms()[1]
                Element_corrected = Poly(
                    [coeff % self.p for coeff in Element.all_coeffs()], x, domain="ZZ"
                )
                Decomposition.append(Element_corrected.as_expr())
            return Decomposition
        else:
            return "Not a generator"


def Affine_encryption(Text_for_encryption: list, Galois: Galois_field, a, b) -> str:
    """
    Encrypts a given text using the affine cipher over the Galois field.

    Args:
        Text_for_encryption: The text to be encrypted
        Galois: The Galois field object
        a: The coefficient of the affine transformation
        b: The coefficient of the affine transformation

    Returns:
        The encrypted text as a string

    """
    a = Poly(a, x)
    b = Poly(b, x)
    Blocks_with_encrypted_text = []
    for i in range(len(Text_for_encryption)):
        # Apply the affine transformation
        Element = a * Text_for_encryption[i] + b
        # Reduce the element modulo the irreducible polynomial
        Element = Poly(Element, x).div(Poly(Galois.irreducible_polynomial, x))[1]
        # Clear the denominator
        Element = Element.clear_denoms()[1]
        # Take the polynomial modulo the prime
        Element_corrected = Poly(
            [coeff % Galois.p for coeff in Element.all_coeffs()], x
        )
        Blocks_with_encrypted_text.append(Element_corrected.as_expr())
    Encrypted_text = "/".join(str(i) for i in Blocks_with_encrypted_text)
    print("The encrypted text as polynomials (Galois field):", Encrypted_text)
    print("The encrypted text as binary:", polynoms_to_binary(Encrypted_text, Galois))
    print(
        "The encrypted text as text:",
        Binary_to_text(polynoms_to_binary(Encrypted_text, Galois)),
    )
    return Encrypted_text


def Affine_decryption(Blocks_with_cipher_text: list, Galois: Galois_field, a, b) -> str:
    """
    Decodes a given text encrypted by the affine cipher over the Galois field.

    Args:
        Blocks_with_cipher_text: A list of polynomials representing the encrypted text
        Galois: A Galois field object
        a: The coefficient of the affine transformation
        b: The coefficient of the affine transformation

    Returns:
        The decoded text as a string
    """
    a = Poly(a, x)
    b = Poly(b, x)

    # Find the modular inverse of the coefficient a
    # Iterate over all elements of the Galois field to find the inverse
    a_inverse = None
    for a_1 in Galois.Galois_field:
        Element = a_1 * a
        Element = Poly(Element, x).div(Poly(Galois.irreducible_polynomial, x))[1]
        Element = Element.clear_denoms()[1]
        Element_corrected = Poly(
            [coeff % Galois.p for coeff in Element.all_coeffs()], x
        ).as_expr()
        if Element_corrected == 1:
            a_inverse = a_1
            break

    # Decrypt each block by applying the inverse affine transformation
    Blocks_with_decryption_text = []
    for block in Blocks_with_cipher_text:
        Element = (sympy.sympify(block) - b) * a_inverse
        Element = Poly(Element, x).div(Poly(Galois.irreducible_polynomial, x))[1]
        Element = Element.clear_denoms()[1]
        Element_corrected = Poly(
            [coeff % Galois.p for coeff in Element.all_coeffs()], x
        )
        Blocks_with_decryption_text.append(Element_corrected.as_expr())

    # Convert decrypted polynomials to binary and then to text
    Open_text_binary = polynoms_to_binary(
        "/".join(str(i) for i in Blocks_with_decryption_text), Galois
    )
    Open_text = Binary_to_text(Open_text_binary)

    # Output the decrypted text
    print("The decrypted text as plain text:", Open_text)
    return Open_text


def polynoms_to_binary(poly_string: str, Galois: Galois_field) -> str:
    """
    Converts a string of polynomials to a binary string.

    Args:
        poly_string: A string of polynomials separated by '/'
        Galois: The Galois field object

    Returns:
        A binary string
    """
    Blocks_with_polynoms = [sympy.sympify(i) for i in poly_string.split("/")]
    Blocks_with_coefficients = get_padded_coeffs(Blocks_with_polynoms, Galois.n)
    Blocks_with_binary = []
    for i in range(len(Blocks_with_coefficients)):
        Blocks_with_binary.append(
            "".join(str(Blocks_with_coefficients[i][j]) for j in range(Galois.n))
        )

    # Join the binary strings into one string
    binary_string = "".join(Blocks_with_binary)
    return binary_string
