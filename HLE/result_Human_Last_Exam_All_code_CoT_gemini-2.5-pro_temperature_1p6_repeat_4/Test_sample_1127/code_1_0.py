import sympy

def solve_minimal_polynomial():
    """
    This function calculates the minimal polynomial of the connective constant mu
    for the given graph G.

    The connective constant mu is derived to be (3 + sqrt(5))/2 using methods
    from statistical mechanics, specifically by calculating the spanning tree
    growth constant of the graph.

    This function starts with this value of mu and finds its minimal polynomial over Q.
    """
    
    # Let 'x' be the variable for the polynomial. In our thoughts, this is mu.
    x = sympy.Symbol('x')

    # The value of the connective constant is mu = (3 + sqrt(5)) / 2.
    # We find the polynomial P(x) = 0 that has mu as a root.
    # x = (3 + sqrt(5)) / 2
    # 2*x = 3 + sqrt(5)
    # 2*x - 3 = sqrt(5)
    
    # Square both sides to eliminate the square root.
    # (2*x - 3)**2 = 5
    # 4*x**2 - 12*x + 9 = 5
    # 4*x**2 - 12*x + 4 = 0
    
    # Divide by 4 to get the monic minimal polynomial.
    # x**2 - 3*x + 1 = 0
    
    # This polynomial is irreducible over the rational numbers because its roots,
    # (3 +/- sqrt(5))/2, are irrational. Thus, it is the minimal polynomial.
    
    polynomial_expr = x**2 - 3*x + 1

    # Extract the coefficients to print them individually as requested.
    poly = sympy.Poly(polynomial_expr, x)
    coeffs = poly.all_coeffs()
    
    c2 = coeffs[0]  # Coefficient of x^2
    c1 = coeffs[1]  # Coefficient of x
    c0 = coeffs[2]  # Constant term
    
    print("The minimal polynomial P(x) for the connective constant is:")
    print(f"P(x) = ({c2})*x^2 + ({c1})*x + ({c0}) = 0")


solve_minimal_polynomial()

# The final answer is the polynomial itself.
print("\n<<<x**2 - 3*x + 1>>>")