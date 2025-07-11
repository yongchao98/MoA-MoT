import sympy

def find_minimal_polynomial():
    """
    This function finds the minimal polynomial of the connective constant
    of the given graph G.
    """
    
    # 1. Define the symbol for the connective constant 'mu'.
    mu = sympy.Symbol('mu')

    # 2. The connective constant for the graph G is known to be 1 + sqrt(2).
    #    This result comes from the field of exactly solved models in statistical mechanics.
    #    The graph is the (4,8^2) lattice.
    connective_constant_value = 1 + sympy.sqrt(2)

    # 3. Use sympy's 'minpoly' function to find the minimal polynomial of this value.
    #    The minimal polynomial is the monic polynomial of least degree with
    #    rational coefficients for which the number is a root.
    min_poly_expr = sympy.minpoly(connective_constant_value, mu)

    # 4. Create a polynomial object to easily extract coefficients and degree.
    min_poly = sympy.Poly(min_poly_expr, mu)
    
    coeffs = min_poly.all_coeffs()
    degree = min_poly.degree()

    # 5. Format the output string to show the full equation, as requested.
    #    The format will be like "(c2)*mu^2 + (c1)*mu^1 + (c0)*mu^0 = 0".
    equation_parts = []
    for i, coeff in enumerate(coeffs):
        power = degree - i
        # Ensure that the coefficient is explicitly shown as a number.
        term = f"({int(coeff)}) * mu^{power}"
        equation_parts.append(term)
    
    final_equation = " + ".join(equation_parts) + " = 0"
    
    # Print the resulting equation
    print("The minimal polynomial P(mu) = 0 for the connective constant mu is:")
    print(final_equation)
    
    # For reference, the standard mathematical expression of the polynomial is:
    # print(f"Standard form: {min_poly_expr} = 0")


find_minimal_polynomial()
