def solve_max_n(k):
    """
    Calculates the maximum value of n for a given k based on the problem's conditions.

    Args:
        k (int): The size of the subsets in the uniform family. Must be an integer > 1.
    """
    if not isinstance(k, int) or k <= 1:
        print("Error: k must be an integer greater than 1.")
        return

    # From the mathematical analysis, the maximum value of n is 2k - 1.
    max_n = 2 * k - 1
    
    print(f"The question is to find the maximum value of n for a given k.")
    print(f"The analysis of the properties of the k-uniform intersecting family shows that the maximum value is obtained when n = 2k - 1.")
    print(f"\nFor the given value k = {k}, the calculation is as follows:")
    
    # We output each number in the final equation as requested.
    print(f"n_max = 2 * {k} - 1 = {max_n}")

# Define the value for k. You can change this to any integer greater than 1.
k_value = 12
solve_max_n(k_value)
