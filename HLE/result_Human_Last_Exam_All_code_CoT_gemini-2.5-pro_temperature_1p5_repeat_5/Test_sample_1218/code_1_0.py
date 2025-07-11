def solve_max_n(k):
    """
    Calculates the maximum value of n in terms of k for a k-uniform
    intersecting family with full differences of size k-1.

    Args:
        k: The size of the subsets in the family (an integer > 1).

    Returns:
        The maximum value of n.
    """
    if not isinstance(k, int) or k <= 1:
        print("Error: k must be an integer greater than 1.")
        return None

    n = 2 * k - 1
    print(f"For k = {k}, the maximum value of n is given by the equation:")
    print(f"n = 2 * {k} - 1 = {n}")
    return n

# Example calculation for k=4
solve_max_n(4)
