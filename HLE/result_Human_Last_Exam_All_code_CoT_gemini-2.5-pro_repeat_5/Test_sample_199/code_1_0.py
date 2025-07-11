import sympy

def solve_minimal_polynomial():
    """
    Calculates and prints the minimal polynomial for the shortest closed geodesic
    on a dodecahedron starting from a vertex.
    """
    # Define the symbolic variable for our polynomial
    x = sympy.Symbol('x')

    # Step 1: Use the known result for the squared length of the shortest path.
    # For a regular dodecahedron with side length 1, the shortest closed
    # geodesic starting from a vertex has a squared length (L^2) of:
    # L^2 = (11 + 3 * sqrt(5)) / 2
    # We are interested in the length L itself.
    l_squared = (11 + 3 * sympy.sqrt(5)) / 2
    length = sympy.sqrt(l_squared)

    # Step 2: Use sympy to find the minimal polynomial of the length L.
    # A minimal polynomial is the simplest polynomial with integer coefficients
    # that has the given number as a root.
    min_poly = sympy.minpoly(length, x)

    # Step 3: Extract the coefficients and print the final equation.
    # The problem asks to show each number in the equation.
    # The polynomial is of the form: a*x^4 + b*x^3 + c*x^2 + d*x + e = 0
    c4 = min_poly.coeff(x, 4)
    c3 = min_poly.coeff(x, 3)
    c2 = min_poly.coeff(x, 2)
    c1 = min_poly.coeff(x, 1)
    c0 = min_poly.coeff(x, 0)

    print("The minimal polynomial for the shortest distance x is P(x) = 0, where:")
    # The f-string formats the final equation, showing each coefficient explicitly.
    print(f"P(x) = {int(c4)}*x^4 + {int(c3)}*x^3 + ({int(c2)})*x^2 + {int(c1)}*x + {int(c0)}")
    print("\nThe equation is:")
    print(f"{int(c4)}*x^4 + {int(c3)}*x^3 + {int(c2)}*x^2 + {int(c1)}*x + {int(c0)} = 0")


solve_minimal_polynomial()