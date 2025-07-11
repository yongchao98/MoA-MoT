def solve_max_n(k):
    """
    Calculates the maximum value of n in terms of k for a k-uniform intersecting family
    with full differences of size k-1.

    Args:
        k (int): The size of the subsets in the family. Must be a positive integer.

    Returns:
        int: The maximum value of n.
    """
    if not isinstance(k, int) or k <= 0:
        raise ValueError("k must be a positive integer.")
    
    # The maximum value of n is given by the formula 2k - 1.
    n = 2 * k - 1
    
    # The problem asks to output the final equation.
    # We use print to display the result as requested.
    print(f"Given k = {k}")
    print(f"The maximum value of n is calculated by the equation: n = 2 * k - 1")
    print(f"n = 2 * {k} - 1")
    print(f"n = {2*k} - 1")
    print(f"n = {n}")
    
    return n

# Example usage:
# Let's use a sample value for k, for example k=3.
k_value = 3
solve_max_n(k_value)

# You can test with other values of k as well
# k_value = 5
# solve_max_n(k_value)