import sympy

def find_closed_form():
    """
    This function finds a closed form for the sum:
    S_n = sum_{k=0 to n} (2k+1)^5 * C(2k,k) * C(2n-2k, n-k)
    
    The method uses generating functions. The sum is the n-th coefficient 
    of the product of two generating functions. The final expression for S_n is 
    4^n * Q_5(n), where Q_5(n) is a polynomial in n. This polynomial is 
    derived and simplified using sympy.
    """

    n = sympy.Symbol('n')

    # This expression for Q_5(n) * 120 is derived from the generating function method.
    # The coefficients 16, 296, 516, 116, 1 are derived from applying the
    # operator (2*x*d/dx + 1)^5 to the generating function of the central
    # binomial coefficients.
    
    poly_times_120 = (
        16 * (n + 1) * n * (n - 1) * (n - 2) * (n - 3)
        + 296 * (n + 2) * (n + 1) * n * (n - 1) * (n - 2)
        + 516 * (n + 3) * (n + 2) * (n + 1) * n * (n - 1)
        + 116 * (n + 4) * (n + 3) * (n + 2) * (n + 1) * n
        + (n + 5) * (n + 4) * (n + 3) * (n + 2) * (n + 1)
    )

    # Simplify the expression to get the polynomial Q_5(n)
    q5_n_poly = sympy.expand(poly_times_120 / 120)

    # Express the resulting polynomial as a fraction with a common denominator
    # for a cleaner representation.
    q5_n_poly_frac = sympy.fraction(sympy.cancel(q5_n_poly))
    numerator = sympy.expand(q5_n_poly_frac[0])
    denominator = q5_n_poly_frac[1]
    
    print("The closed form for the sum is:")
    
    # We will output the formula in the form S_n = 4**n * (polynomial expression)
    # The user is asked to have each number in the final equation. We achieve this
    # by printing the expanded polynomial coefficients.

    num_poly = sympy.Poly(numerator, n)
    num_coeffs = num_poly.all_coeffs()
    num_terms = []
    for i, coeff in enumerate(num_coeffs):
        power = len(num_coeffs) - 1 - i
        if coeff == 0:
            continue
        term = f"{coeff}*n**{power}"
        num_terms.append(term)
    
    numerator_str = " + ".join(num_terms).replace("+ -", "- ")
    
    print(f"S_n = 4**n * ( ( {numerator_str} ) / {denominator} )")

find_closed_form()