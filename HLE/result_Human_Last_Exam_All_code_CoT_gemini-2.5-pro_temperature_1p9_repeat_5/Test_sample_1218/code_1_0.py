def solve_max_n(k):
    """
    Calculates the maximum value of n for a given k based on the derived formula.
    
    The problem is to find the maximum value of n such that there exists a
    k-uniform intersecting family F on a set of n elements with full differences
    of size k-1. The derivation shows this value is n = k^2 - k + 1.

    Args:
        k (int): The size of the subsets in the uniform family (must be >= 2).
    """
    if not isinstance(k, int) or k < 2:
        print("Error: k must be an integer greater than or equal to 2.")
        return

    # The formula for the maximum value of n is k^2 - k + 1.
    n_max = k**2 - k + 1

    print(f"For k = {k}, the maximum value of n is determined by the formula n = k^2 - k + 1.")
    print("The calculation is as follows:")
    
    # As requested, output each number in the final equation.
    print(f"n = ({k})^2 - ({k}) + 1 = {k*k} - {k} + 1 = {n_max}")
    print(f"\nFinal Answer: The maximum value of n is {n_max}.")

# Example: solve the problem for a sample value of k=5
solve_max_n(5)