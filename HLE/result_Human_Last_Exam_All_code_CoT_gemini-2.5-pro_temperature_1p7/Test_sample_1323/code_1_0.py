import sympy

def solve_expression():
    """
    This function represents the derived expression for '?_1' using sympy
    and prints it to the console.
    """

    # Define the symbolic variables and functions
    # x is a vector (x1, x2), and h is a function of x.
    x1, x2 = sympy.symbols('x1 x2')
    h = sympy.Function('h')(x1, x2)

    # i and j are the indices for the partial derivatives
    i, j = sympy.symbols('i j', integer=True)

    # The result involves the Kronecker Delta function
    delta_ij = sympy.KroneckerDelta(i, j)

    # The derived numerical constant is 1/2
    # The final equation is ?_1 = (1/2) * h(x) * delta_ij
    # We output the number 1/2 here as part of the equation construction
    constant = sympy.Rational(1, 2)

    # Construct the final expression for ?_1
    q1_expression = constant * h * delta_ij

    # Print the general result
    print("The derived expression for ?_1 is:")
    sympy.pprint(q1_expression, use_unicode=True)
    print("\nwhere δ(i, j) is the Kronecker delta.")
    
    # Show the specific expressions for the components of the 2x2 matrix
    print("\nIn matrix form, where i, j are row and column indices (1 or 2):")
    matrix_expr = sympy.Matrix([[q1_expression.subs({i:r, j:c}) for c in [1,2]] for r in [1,2]])
    sympy.pprint(matrix_expr, use_unicode=True)


if __name__ == '__main__':
    solve_expression()