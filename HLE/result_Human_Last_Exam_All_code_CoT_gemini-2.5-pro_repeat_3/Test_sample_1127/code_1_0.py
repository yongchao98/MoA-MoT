import sympy

def solve():
    """
    This function calculates and prints the minimal polynomial for the
    connective constant of the specified graph.
    """
    
    # The connective constant 'mu' for the given graph (a "brick-wall" lattice)
    # is known to be sqrt(3). We represent this using sympy.
    mu_val = sympy.sqrt(3)
    
    # We define a symbol 'x' to represent the variable in our polynomial.
    x = sympy.Symbol('x')
    
    # We use sympy's 'minpoly' function to find the minimal polynomial of mu_val.
    # This function finds the monic polynomial with rational coefficients of the
    # smallest degree that has mu_val as a root.
    min_poly = sympy.minpoly(mu_val, x)
    
    # The minimal polynomial is x**2 - 3. We will now print the equation P(mu) = 0.
    
    # Get the polynomial as a Poly object to easily access coefficients.
    poly_obj = min_poly.as_poly(x)
    
    # Get all coefficients from the highest degree to the constant term.
    # For x**2 - 3, this will be [1, 0, -3].
    coeffs = poly_obj.all_coeffs()
    degree = poly_obj.degree()
    
    # Build the equation string, showing each number as requested.
    equation_parts = []
    for i, coeff in enumerate(coeffs):
        power = degree - i
        
        # Add sign for terms after the first one
        if i > 0:
            if coeff < 0:
                equation_parts.append(" - ")
                # Use the absolute value for printing
                coeff = -coeff
            else:
                equation_parts.append(" + ")
        
        # Append the coefficient
        equation_parts.append(str(coeff))
        
        # Append the variable part " * mu^power "
        if power > 1:
            equation_parts.append(f" * mu^{power}")
        elif power == 1:
            equation_parts.append(" * mu")
        # No variable part for the constant term (power = 0)
        
    equation_parts.append(" = 0")
    
    final_equation = "".join(equation_parts)
    print(final_equation)

solve()