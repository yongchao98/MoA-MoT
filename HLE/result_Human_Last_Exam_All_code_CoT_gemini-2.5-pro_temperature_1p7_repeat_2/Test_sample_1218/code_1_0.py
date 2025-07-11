def solve_max_n(k: int):
    """
    Calculates the maximum value of n for a k-uniform intersecting family
    with full differences of size k-1.

    Args:
        k: The size of the subsets in the family F (must be an integer >= 2).
    """
    if not isinstance(k, int) or k < 2:
        print("Error: k must be an integer greater than or equal to 2.")
        return

    # Based on the analysis, the maximum value of n is 2k - 1.
    n = 2 * k - 1

    print(f"For a given k-uniform intersecting family with full differences of size k-1,")
    print(f"the maximum value of n in terms of k is given by the equation:")
    print(f"n_max = 2 * k - 1")
    print()
    print(f"For k = {k}, the calculation is:")
    # The prompt requests to "output each number in the final equation",
    # which is interpreted as showing the steps of the calculation.
    print(f"n_max = 2 * {k} - 1")
    print(f"n_max = {2 * k} - 1")
    print(f"n_max = {n}")


# As an example, let's use a value for k.
# You can change this value to test other cases.
k_example = 5
solve_max_n(k_example)
