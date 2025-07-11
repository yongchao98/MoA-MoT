import math

def find_subset():
    """
    Finds the subset of integers t from the given set for which any n x n grid
    has an even number of t-omino tilings.

    The property holds if and only if t is not a perfect square.
    """
    t_values = [2, 3, 4, 5, 7, 9, 15]
    result_subset = []

    for t in t_values:
        # Check if t is a perfect square.
        # A number is a perfect square if its square root is an integer.
        sqrt_t = math.isqrt(t)
        if sqrt_t * sqrt_t != t:
            result_subset.append(t)

    print("The subset of integers for which the statement is true is:")
    print(result_subset)

find_subset()