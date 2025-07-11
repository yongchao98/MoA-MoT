def solve_max_n(k: int):
    """
    Calculates the maximum value of n for a k-uniform intersecting family
    with full differences of size k-1.

    Args:
        k: The size of the subsets in the family.
    """
    if k <= 1:
        print("The value of k must be an integer greater than 1.")
        return

    # Based on combinatorial analysis, the maximum value of n is 2k - 1.
    # We will demonstrate this for the given value of k.
    
    n_max = 2 * k - 1
    
    # Print the equation as requested
    print(f"For k = {k}, the maximum value of n is given by the equation:")
    print(f"n_max = 2 * {k} - 1 = {n_max}")


# --- Main execution ---
# We will use a sample value for k, for example, k=10.
k_value = 10
solve_max_n(k_value)