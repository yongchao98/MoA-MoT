import sympy

def find_closed_form():
    """
    This function calculates the closed form for the given sum using symbolic mathematics.
    The derivation involves generating functions, leading to an expression for the sum
    as 4^n times a polynomial in n. This script computes the coefficients of that polynomial.
    """
    n = sympy.symbols('n')

    # From the generating function derivation, the sum S_n can be expressed as:
    # S_n = 4^n * Q_5(n), where Q_5(n) is a polynomial in n.
    # Q_5(n) is given by the following sum of binomial coefficients:
    # This comes from S_n = [x^n] P_5(x) / (1-4x)^6, where P_5(x) is a known polynomial.
    # P_5(x) = 4096x^4 + 18944x^3 + 8256x^2 + 464x + 1
    # Q_5(n) = sum_{j=0 to 4} c_j * 4^(-j) * binomial(n-j+5, 5)
    # where c_j are coefficients of P_5(x)
    
    Q5_n = (1 * sympy.binomial(n + 5, 5) * sympy.Rational(1, 4**0) +
            464 * sympy.binomial(n + 4, 5) * sympy.Rational(1, 4**1) +
            8256 * sympy.binomial(n + 3, 5) * sympy.Rational(1, 4**2) +
            18944 * sympy.binomial(n + 2, 5) * sympy.Rational(1, 4**3) +
            4096 * sympy.binomial(n + 1, 5) * sympy.Rational(1, 4**4))

    # Expand the binomial coefficients and simplify the expression to get the polynomial Q_5(n)
    poly_Q5_n = sympy.expand(Q5_n)
    
    # The polynomial is a rational function. Get its numerator and denominator.
    num, den = sympy.fraction(poly_Q5_n)
    
    # The numerator is the polynomial part of our answer.
    num_poly = sympy.Poly(num, n)
    
    # The final closed form is 4^n * (Numerator / Denominator)
    # We will print the coefficients of the numerator polynomial and the denominator.
    
    print("The closed form for the sum is:")
    print("S_n = 4^n * ( (A*n^5 + B*n^4 + C*n^3 + D*n^2 + E*n) / G )")
    print("\nwhere the coefficients are:")
    print(f"A = {num_poly.coeff_monomial(n**5)}")
    print(f"B = {num_poly.coeff_monomial(n**4)}")
    print(f"C = {num_poly.coeff_monomial(n**3)}")
    print(f"D = {num_poly.coeff_monomial(n**2)}")
    print(f"E = {num_poly.coeff_monomial(n**1)}")
    print(f"G = {den}")

if __name__ == '__main__':
    find_closed_form()