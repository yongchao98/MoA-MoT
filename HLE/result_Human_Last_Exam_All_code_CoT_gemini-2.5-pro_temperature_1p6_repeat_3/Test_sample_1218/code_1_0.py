def solve_max_n(k):
    """
    Calculates the maximum value of n in terms of k for a k-uniform intersecting family
    with full differences of size k-1.

    The maximum value is given by the formula n = k^2 - k + 1.

    Args:
        k (int): The size of the subsets in the family (k-uniform). Must be >= 2.
    """
    if not isinstance(k, int) or k < 2:
        print("Error: k must be an integer greater than or equal to 2.")
        return

    # According to extremal set theory, the maximum value of n is k^2 - k + 1.
    # Step 1: Calculate k squared
    k_squared = k * k
    
    # Step 2: Calculate k squared minus k
    k_squared_minus_k = k_squared - k
    
    # Step 3: Add 1 to get the final result for n
    n = k_squared_minus_k + 1
    
    # Print the equation step-by-step
    print(f"For k = {k}, the maximum value of n is calculated as follows:")
    print(f"n = k^2 - k + 1")
    print(f"n = {k}^2 - {k} + 1")
    print(f"n = {k_squared} - {k} + 1")
    print(f"n = {k_squared_minus_k} + 1")
    print(f"n = {n}")

# Example usage:
# Let's use an integer k as an example to demonstrate the calculation.
# For instance, if k=3:
k_example = 3
solve_max_n(k_example)

# Let's try another example, k=4:
# k_example_2 = 4
# solve_max_n(k_example_2)
