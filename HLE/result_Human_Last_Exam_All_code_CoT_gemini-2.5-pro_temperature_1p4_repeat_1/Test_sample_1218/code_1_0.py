def solve_max_n(k):
    """
    Calculates the maximum value of n in terms of k for the given problem.

    The derivation shows that n must be at least 2k - 1,
    and that n cannot be 2k or greater.
    Thus, the maximum value for n is 2k - 1.

    Args:
        k (int): The size of the subsets in the uniform family. Should be an integer >= 2.
    """
    if not isinstance(k, int) or k < 2:
        print("Please provide an integer k where k >= 2.")
        return

    n_max = 2 * k - 1

    print(f"For k = {k}, the maximum value of n is determined by the equation:")
    # The prompt requests to output each number in the final equation.
    print(f"n_max = 2 * {k} - 1")
    print(f"n_max = {2 * k} - 1")
    print(f"n_max = {n_max}")


# Example usage:
# You can change this value to test with other k.
k_value = 4
solve_max_n(k_value)
