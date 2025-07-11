import math

def calculate_max_n_for_intersecting_family(k):
    """
    This function calculates the maximum value of n for a k-uniform
    intersecting family with full differences of size k-1.

    Args:
        k (int): The size of subsets, an integer >= 2.
    """
    if not isinstance(k, int) or k < 2:
        print("Error: k must be an integer greater than or equal to 2.")
        return

    # Based on the derivation, the maximum value of n is 2k - 1.
    n = 2 * k - 1

    print(f"For k = {k}, the maximum possible value of n is calculated as follows:")
    
    # Printing the final equation with each number explicitly.
    print("\nEquation:")
    print(f"  n = 2 * k - 1")
    
    print("\nSubstitution:")
    print(f"  n = {2} * {k} - {1}")
    
    print("\nResult:")
    print(f"  n = {n}")


# Example: Calculate the maximum n for k=5
k_value = 5
calculate_max_n_for_intersecting_family(k_value)
