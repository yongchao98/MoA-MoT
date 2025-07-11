def solve_max_n(k):
    """
    Calculates the maximum value of n in terms of k for the given combinatorial problem.

    Args:
        k (int): The size of the subsets in the family, must be an integer >= 2.

    Returns:
        None. Prints the result.
    """
    if not isinstance(k, int) or k < 2:
        print("Error: k must be an integer greater than or equal to 2.")
        return

    # The maximum value of n is given by the formula n = 2k - 1.
    # We established this by:
    # 1. Proving that n must be at least 2k - 1 for the difference sets to be formed.
    # 2. Showing that n = 2k - 1 works with the family F = all k-subsets of [2k-1].
    # 3. Proving that for n >= 2k, any such family would be forced to contain two
    #    disjoint sets, which contradicts the "intersecting" property.

    n = 2 * k - 1

    print(f"For a k-uniform intersecting family with full differences of size k-1:")
    print(f"Given k = {k}")
    print("The maximum value of n is calculated using the formula: n = 2 * k - 1")
    print(f"n = 2 * {k} - 1")
    print(f"n = {2*k} - 1")
    print(f"n = {n}")

# Example usage with a specific value for k.
# You can change this value to see the result for a different k.
k_value = 4
solve_max_n(k_value)