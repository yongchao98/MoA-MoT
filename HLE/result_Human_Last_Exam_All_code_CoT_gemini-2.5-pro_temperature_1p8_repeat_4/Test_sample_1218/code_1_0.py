def solve_max_n(k):
    """
    Calculates the maximum value of n for a k-uniform intersecting family
    with full differences of size k-1.

    Args:
        k (int): The size of the subsets in the family. k must be an integer >= 2.

    Returns:
        int: The maximum value of n.
    """
    if not isinstance(k, int) or k < 2:
        print("Error: k must be an integer greater than or equal to 2.")
        return None

    # The maximum value of n is given by the formula: n = 2k - 1
    max_n = 2 * k - 1

    print(f"For k = {k}, the maximum value of n is given by the equation:")
    print(f"n_max = 2 * k - 1")
    print(f"n_max = 2 * {k} - 1 = {max_n}")
    return max_n

# Example usage with k=4
k_value = 4
solve_max_n(k_value)

print("-" * 20)

# Example usage with k=10
k_value_2 = 10
solve_max_n(k_value_2)
