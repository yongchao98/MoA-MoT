def solve_max_n(k):
    """
    Calculates the maximum value of n for a given k based on the formula n = k^2 - k + 1.

    Args:
        k (int): The size of the subsets in the family F. Must be >= 2.
    """
    if not isinstance(k, int) or k < 2:
        print("Error: k must be an integer greater than or equal to 2.")
        return

    # The maximum value of n is given by the formula k^2 - k + 1.
    n_max = k**2 - k + 1
    k_squared = k**2

    print(f"For k = {k}, the maximum value of n is given by the formula: n = k^2 - k + 1")
    print("The calculation is as follows:")
    print(f"n = {k}^2 - {k} + 1")
    print(f"n = {k_squared} - {k} + 1")
    print(f"n = {n_max}")

# Example usage of the function with k=3.
# For k=3, a projective plane of order 2 (the Fano plane) exists.
# n = 3^2 - 3 + 1 = 7.
solve_max_n(3)