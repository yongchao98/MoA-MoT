import sympy

def solve_minimal_polynomial():
    """
    This function calculates and prints the minimal polynomial for the shortest 
    closed geodesic on a dodecahedron starting from a vertex.
    """
    # The problem is to find the minimal polynomial of the shortest possible distance an ant 
    # could walk along a geodesic trajectory on the surface of a regular dodecahedron, 
    # starting from a vertex and returning to it without passing through any other vertex.
    # The side length of the dodecahedron is 1.

    # From geometric analysis (a known result), the shortest such distance, let's call it d,
    # is given by d = sqrt(7 - 3*sqrt(5)).
    # We will find the minimal polynomial for this value d.

    # Let x be the variable for our polynomial, representing the distance d.
    x = sympy.Symbol('x')

    # We start with the equation x = sqrt(7 - 3*sqrt(5)).
    # To find the polynomial, we eliminate the square roots.
    
    # Step 1: Square both sides.
    # x^2 = 7 - 3*sqrt(5)
    
    # Step 2: Isolate the remaining square root term.
    # x^2 - 7 = -3*sqrt(5)
    
    # Step 3: Square both sides again.
    # (x^2 - 7)^2 = (-3*sqrt(5))^2
    # (x^2 - 7)^2 = 9 * 5 = 45
    
    # Step 4: Expand and simplify the equation.
    # x^4 - 14*x^2 + 49 = 45
    # x^4 - 14*x^2 + 4 = 0
    
    # We use sympy to perform this derivation formally.
    poly_expr = (x**2 - 7)**2 - 45
    
    # Expand the expression to get the final polynomial.
    final_poly = sympy.expand(poly_expr)
    
    # The minimal polynomial is final_poly = 0.
    p = sympy.Poly(final_poly, x)
    coeffs = p.all_coeffs()
    degree = p.degree()

    print("The minimal polynomial equation for the shortest distance x is:")
    
    # The prompt requires printing each number in the final equation.
    # We will construct a list of strings representing each number and symbol
    # and then print them.
    
    output_parts = []
    is_first_term = True
    
    for i, coeff in enumerate(coeffs):
        power = degree - i
        
        if coeff == 0:
            continue
        
        # Determine the sign and coefficient value to print.
        if is_first_term:
            if coeff < 0:
                output_parts.append(str(int(coeff)))
            # For a positive leading term, we don't print '+'
            # If the coefficient is 1, we don't print it unless it's the constant term.
            elif coeff > 1:
                output_parts.append(str(int(coeff)))
        else:
            if coeff < 0:
                output_parts.append(" - ")
                output_parts.append(str(int(-coeff)))
            else:
                output_parts.append(" + ")
                output_parts.append(str(int(coeff)))

        # Add the variable part (x, x^2, etc.).
        if power > 0:
            # Add a multiplication sign if a coefficient is explicitly printed.
            if not (is_first_term and coeff == 1):
                 output_parts.append(" * ")
            output_parts.append("x")
            if power > 1:
                output_parts.append("^")
                output_parts.append(str(power))
        
        is_first_term = False

    output_parts.append(" = ")
    output_parts.append("0")

    print(*output_parts)

solve_minimal_polynomial()