def solve_for_n(k):
    """
    This function calculates the maximum value of n in terms of k
    for which a k-uniform intersecting family with full differences
    of size k-1 exists.

    The formula derived is n = 2k - 1.
    """
    if k < 1:
        print("k must be a positive integer.")
        return

    # Calculate n using the formula
    n = 2 * k - 1

    # Output the result, showing each number in the equation
    print(f"For k = {k}, the maximum value of n is given by the equation:")
    print(f"n = 2 * {k} - 1")
    print(f"n = {2 * k} - 1")
    print(f"n = {n}")

# Example usage with k = 10
k_example = 10
solve_for_n(k_example)