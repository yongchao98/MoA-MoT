def solve_max_n(k):
    """
    Calculates the maximum value of n for a given k based on the formula n = 2k - 1.
    Prints the formula and the result for the given k.
    
    Args:
        k (int): The size of the subsets in the family F, must be an integer >= 2.
    """
    if not isinstance(k, int) or k < 2:
        print("Error: k must be an integer greater than or equal to 2.")
        return

    # The maximum value of n is given by the formula n = 2k - 1.
    n = 2 * k - 1
    
    print(f"For a k-uniform intersecting family with full differences of size k-1, the maximum value of n is given by the formula:")
    print("n = 2*k - 1")
    print("\nFor the specific case where k is an integer, for example k =", k, ", the calculation is:")
    print("n = 2 *", k, "- 1 =", n)

# Example execution with a value for k.
# You can change this value to see the result for a different k.
example_k = 4
solve_max_n(example_k)