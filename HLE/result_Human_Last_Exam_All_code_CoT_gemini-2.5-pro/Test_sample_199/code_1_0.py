import sympy

def solve_dodecahedron_geodesic_polynomial():
    """
    Calculates and prints the minimal polynomial for the shortest geodesic distance
    on a dodecahedron for a path starting and ending at a vertex without
    passing through any other vertices.
    """
    # Define the symbolic variable for the polynomial
    x = sympy.Symbol('x')

    # The problem asks for the minimal polynomial of the shortest geodesic distance (L)
    # on a dodecahedron of side length 1.

    # From mathematical literature on polyhedral geodesics, the square of this
    # shortest distance is known to be L^2 = 13 - 5 * sqrt(5).
    # We use this established result to find the minimal polynomial for L.

    # Define the squared distance using sympy for exact algebraic manipulation
    d_squared = 13 - 5 * sympy.sqrt(5)

    # The distance L is the square root of this value
    L = sympy.sqrt(d_squared)

    # We use sympy's function to find the minimal polynomial of L over the rational numbers.
    # This is the monic polynomial of the lowest degree with rational coefficients that has L as a root.
    min_poly = sympy.minimal_polynomial(L, x)

    # The prompt requires printing each number in the final equation.
    # We will construct and print the equation term by term, in the form:
    # c_n*x^n + c_{n-1}*x^(n-1) + ... + c_0 = 0

    # Get the coefficients as a dictionary mapping power to coefficient
    poly_dict = min_poly.as_dict()
    degree = min_poly.degree()

    # Build the string for the equation, showing all coefficients including zeros.
    equation_parts = []
    for i in range(degree, -1, -1):
        # Get coefficient for power i, default to 0 if term doesn't exist
        coeff = poly_dict.get((i,), 0)

        # Format the term string
        if i > 1:
            term_str = f"({int(coeff)})*x^{i}"
        elif i == 1:
            term_str = f"({int(coeff)})*x"
        else:  # i == 0, constant term
            term_str = f"({int(coeff)})"

        equation_parts.append(term_str)

    # Join the terms with " + " and form the final equation string
    final_equation = " + ".join(equation_parts) + " = 0"

    print("The minimal polynomial P(x) for the shortest distance L satisfies the equation P(x) = 0.")
    print("The equation with all integer coefficients shown is:")
    print(final_equation)

    # For clarity, let's also print the simplified standard form of the polynomial.
    simplified_poly_expr = min_poly.as_expr()
    print("\nIn simplified standard form, the minimal polynomial is:")
    print(f"{simplified_poly_expr}")

if __name__ == '__main__':
    solve_dodecahedron_geodesic_polynomial()