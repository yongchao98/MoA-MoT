def solve_max_n(k):
    """
    Calculates the maximum value of n for a k-uniform intersecting family
    with full differences of size k-1.

    The formula is derived from combinatorial arguments in extremal set theory.
    The maximum value is n = k^2 - k + 1.

    Args:
        k (int): The size of the subsets in the family. Must be >= 2.

    Returns:
        None. Prints the result.
    """
    if not isinstance(k, int) or k < 2:
        print("Error: k must be an integer greater than or equal to 2.")
        return

    # Calculate n using the formula n = k^2 - k + 1
    n = k**2 - k + 1

    # Print the explanation and the final equation
    print(f"For a k-uniform intersecting family with full differences of size k-1,")
    print(f"the maximum value of n in terms of k is given by the formula: n = k^2 - k + 1")
    print("\nFor k = {}, the calculation is:".format(k))
    print(f"n = {k}^2 - {k} + 1")
    print(f"n = {k**2} - {k} + 1")
    print(f"n = {n}")


# Example usage with k = 3
k_value = 3
solve_max_n(k_value)

# Example usage with k = 4
# print("\n" + "="*20 + "\n")
# k_value = 4
# solve_max_n(k_value)