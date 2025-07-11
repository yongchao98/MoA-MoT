import sympy

def solve_minimal_polynomial():
    """
    This function calculates and prints the minimal polynomial for the shortest
    closed geodesic on a dodecahedron starting and ending at a vertex.
    """
    
    # Define x as a symbolic variable for the polynomial
    x = sympy.symbols('x')

    # The shortest distance, d, is sqrt(5 + 2*sqrt(5)) for a dodecahedron of side length 1.
    # We create this expression using sympy.
    d_expression = sympy.sqrt(5 + 2 * sympy.sqrt(5))

    print("The problem is to find the minimal polynomial of the shortest geodesic distance 'd' on a dodecahedron.")
    print(f"The value of this distance is d = sqrt(5 + 2*sqrt(5)).")
    print(f"Numerically, d is approximately: {d_expression.evalf()}")
    print("-" * 30)

    # Sympy's minimal_polynomial function can find the polynomial directly.
    min_poly = sympy.minimal_polynomial(d_expression, x)

    print("The minimal polynomial P(x) for the distance d is:")
    sympy.pprint(min_poly)
    print("-" * 30)
    
    # Extract the coefficients to display the final equation as requested.
    # The polynomial is of the form c4*x^4 + c2*x^2 + c0 = 0
    c4 = min_poly.coeff(x, 4)
    c3 = min_poly.coeff(x, 3)
    c2 = min_poly.coeff(x, 2)
    c1 = min_poly.coeff(x, 1)
    c0 = min_poly.coeff(x, 0)
    
    print("The final equation P(d) = 0, for which d is a root, is structured as follows:")
    # We explicitly show each coefficient, even if it's zero.
    print(f"({c4})*d^4 + ({c3})*d^3 + ({c2})*d^2 + ({c1})*d + ({c0}) = 0")

solve_minimal_polynomial()
