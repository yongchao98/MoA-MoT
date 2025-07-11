def solve_max_n(k):
    """
    Calculates the maximum value of n for a given k based on the derivation.
    The formula is n = 2*k - 1.
    """
    if not isinstance(k, int) or k < 1:
        print("Error: k must be a positive integer.")
        return

    # Calculate n using the formula
    n = 2 * k - 1
    
    # Print the explanation and the final equation as requested
    print(f"For k = {k}, the maximum value of n such that there exists a {k}-uniform intersecting family")
    print("with full differences of size k-1 is given by the formula n = 2*k - 1.")
    print("\nCalculation:")
    print(f"n = 2 * {k} - 1 = {n}")

# Example usage with k=4
k_example = 4
solve_max_n(k_example)