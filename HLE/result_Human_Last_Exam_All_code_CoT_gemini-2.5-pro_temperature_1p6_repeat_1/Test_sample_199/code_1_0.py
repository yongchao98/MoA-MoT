import sympy

def solve_ant_path():
    """
    Calculates the minimal polynomial for the shortest closed geodesic on a
    dodecahedron starting and ending at a vertex.
    """
    # Step 1: Define the expression for the shortest distance squared.
    # This value, d^2 = 5 + 2*sqrt(5), is a known result from the geometry of 
    # geodesics on a dodecahedron with side length 1. The path is a closed
    # loop starting and ending at a vertex without passing through any others.
    d_squared = 5 + 2 * sympy.sqrt(5)

    # Step 2: The distance 'd' is the square root of this expression.
    d = sympy.sqrt(d_squared)

    # Step 3: Use sympy to find the minimal polynomial of 'd' with 'x' as the variable.
    # A minimal polynomial is the polynomial of least degree with rational 
    # coefficients and a leading coefficient of 1 that has 'd' as a root.
    x = sympy.Symbol('x')
    min_poly = sympy.minimal_polynomial(d, x)

    # Step 4: Display the final equation and its coefficients.
    # The minimal polynomial found is x**4 - 10*x**2 + 5.
    # We will format this into a readable equation P(x) = 0.
    
    c4 = int(min_poly.coeff(x, 4))
    c3 = int(min_poly.coeff(x, 3))
    c2 = int(min_poly.coeff(x, 2))
    c1 = int(min_poly.coeff(x, 1))
    c0 = int(min_poly.coeff(x, 0))

    print("The minimal polynomial equation for the distance 'x' is:")
    print(f"{c4}*x^4 + ({c2})*x^2 + {c0} = 0")

    print("\nThe numbers (coefficients) in the final equation are:")
    print(f"Coefficient of x^4: {c4}")
    print(f"Coefficient of x^3: {c3}")
    print(f"Coefficient of x^2: {c2}")
    print(f"Coefficient of x^1: {c1}")
    print(f"Constant term (c0): {c0}")

solve_ant_path()