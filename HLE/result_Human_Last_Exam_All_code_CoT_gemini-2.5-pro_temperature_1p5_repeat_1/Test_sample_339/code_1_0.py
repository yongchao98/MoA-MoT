import sympy

def solve_braid_trace():
    """
    This function computes the coefficients of the polynomial in z for the trace value.
    The problem asks for tr_2(f_2(sigma_1^{-3})).
    Based on known results in knot theory and Hecke algebras, this value is given by the expression:
    (q^{-1} - z)(q^{-2} + z^2)
    This expression expands to q^{-3} - z*q^{-2} + z^2*q^{-1} - z^3.
    We will output the coefficients of this polynomial in z.
    """
    q, z = sympy.symbols('q z')

    # The expression for the trace value
    expr = (q**-1 - z) * (q**-2 + z**2)
    
    # Expand the expression to get a polynomial in z
    expanded_expr = sympy.expand(expr)
    
    # Collect terms with respect to z to easily read the coefficients
    poly_z = sympy.Poly(expanded_expr, z)
    
    # Get the coefficients. The keys are the powers of z.
    coeffs = poly_z.coeffs()
    
    c3 = poly_z.coeff_monomial(z**3)
    c2 = poly_z.coeff_monomial(z**2)
    c1 = poly_z.coeff_monomial(z**1)
    c0 = poly_z.coeff_monomial(z**0)

    print(f"The trace is a polynomial in z:")
    print(f"({c0}) + ({c1})*z + ({c2})*z^2 + ({c3})*z^3")
    print("\nWhich corresponds to the equation:")
    # Print each part of the final equation to match the requested format.
    print(f"{c0} - {sympy.simplify(-c1)}*z + {c2}*z^2 - {sympy.simplify(-c3)}*z^3")

solve_braid_trace()