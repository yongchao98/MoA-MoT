import sympy

def solve_ant_path():
    """
    This script calculates the minimal polynomial of the shortest geodesic distance
    an ant can walk on a regular dodecahedron (side length 1) starting from a vertex 
    and returning to it without passing through any other vertex.
    """

    # Step 1: Define the squared shortest distance in terms of algebraic numbers.
    # The shortest path corresponds to a (1,1) vector on the unfolded surface net.
    # Its squared length d^2 is given by d^2 = 2 + 2*cos(108 degrees).
    # We use the golden ratio phi for an exact representation.
    # phi = (1 + sqrt(5))/2
    # cos(108) = (1 - phi)/2
    # So, d^2 = 2 + 2 * (1 - phi)/2 = 3 - phi.

    # We use sympy to handle the algebraic numbers.
    phi = (1 + sympy.sqrt(5)) / 2
    d_squared = 3 - phi

    # Step 2: Define the distance 'd' and find its minimal polynomial.
    # Let x be the variable for the polynomial.
    x = sympy.Symbol('x')
    # The distance d is the square root of d_squared.
    d = sympy.sqrt(d_squared)

    # The minimal_polynomial function finds the polynomial with rational coefficients
    # of the lowest degree that has 'd' as a root.
    min_poly = sympy.minimal_polynomial(d, x)

    # Step 3: Output the result as requested.
    # The problem asks to output each number in the final equation.
    # We will print the polynomial equation in a readable format and then list its coefficients.

    coeffs = min_poly.all_coeffs()
    degree = min_poly.degree()

    equation_str = ""
    for i, coeff in enumerate(coeffs):
        power = degree - i
        
        # Skip terms with zero coefficient
        if coeff == 0:
            continue

        # Determine the sign for the term
        if i == 0:
            # First term sign logic
            if coeff < 0:
                sign = "- "
            else:
                sign = ""
        else:
            # Subsequent terms sign logic
            if coeff < 0:
                sign = " - "
            else:
                sign = " + "
        
        # Use the absolute value of the coefficient for printing
        abs_coeff = abs(coeff)

        # Format the coefficient string (don't print 1 for non-constant terms)
        if abs_coeff == 1 and power > 0:
            coeff_str = ""
        else:
            coeff_str = str(abs_coeff)

        # Format the variable part of the term
        if power == 1:
            var_str = "x"
        elif power == 0:
            var_str = ""
        else:
            var_str = f"x^{power}"

        # Add a space between coefficient and variable if both are present
        if coeff_str and var_str:
            term_space = " "
        else:
            term_space = ""

        equation_str += f"{sign}{coeff_str}{term_space}{var_str}"

    print("The minimal polynomial equation for the shortest distance x is:")
    print(f"{equation_str} = 0")
    print("\nThe numbers (coefficients) of the polynomial P(x) = a*x^4 + b*x^3 + c*x^2 + d*x + e are:")
    for c in coeffs:
        print(int(c))

if __name__ == '__main__':
    solve_ant_path()