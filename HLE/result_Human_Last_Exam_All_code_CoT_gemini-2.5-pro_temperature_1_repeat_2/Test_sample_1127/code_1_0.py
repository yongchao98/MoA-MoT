def solve_connective_constant_polynomial():
    """
    This function provides the minimal polynomial for the connective constant
    of the specified graph G.

    The derivation is based on a known result from statistical mechanics.
    The connective constant, mu, is sqrt(4 + 2*sqrt(7)).
    From this, we derive the polynomial P(x) = 0 that mu satisfies.
    P(x) = x^4 - 8*x^2 - 12.
    """

    # Define the coefficients of the minimal polynomial P(x) = 0.
    # The polynomial is P(x) = c4*x^4 + c3*x^3 + c2*x^2 + c1*x + c0.
    c4 = 1
    c3 = 0
    c2 = -8
    c1 = 0
    c0 = -12

    # The minimal polynomial equation is x^4 - 8*x^2 - 12 = 0.
    # The prompt asks to output each number in the final equation.
    # We will print the equation and then list the non-zero coefficients,
    # which are the numbers in the simplified final equation.

    equation_str = "x^4 - 8*x^2 - 12 = 0"

    print("The minimal polynomial (over Q) of the connective constant is P(x) = 0, where:")
    print(equation_str)

    print("\nThe numbers (non-zero coefficients) in this equation are:")
    print(f"The coefficient of x^4 is: {c4}")
    print(f"The coefficient of x^2 is: {c2}")
    print(f"The constant term is: {c0}")

solve_connective_constant_polynomial()