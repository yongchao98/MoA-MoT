import numpy as np

def solve_for_lambda(n):
    """
    Calculates the values of λ for which the given integral equation has no solutions.

    Args:
        n (int): The integer exponent in the integral equation. n must be >= 1.

    Returns:
        numpy.ndarray: An array of real values for λ.
    """
    if not isinstance(n, int) or n < 1:
        raise ValueError("n must be an integer greater than or equal to 1.")

    # Construct the n x n matrix A.
    # A_mk = 1 / (m + n - k) where m, k are 0-indexed.
    A = np.zeros((n, n))
    for m in range(n):
        for k in range(n):
            A[m, k] = 1.0 / (m + n - k)

    # Find the eigenvalues of A.
    # The operator is self-adjoint, so its eigenvalues are real.
    # The matrix A is a representation in a non-orthogonal basis, but its eigenvalues
    # are the same as the operator's eigenvalues and are therefore real.
    eigenvalues = np.linalg.eigvals(A)
    
    # The values of λ for which no solution exists are the reciprocals of the eigenvalues.
    # We filter out potential zero eigenvalues which would lead to infinite lambda.
    # However, for this problem, the matrix A is known to be invertible.
    lambda_values = 1.0 / eigenvalues

    return np.sort(lambda_values)

# The user did not specify a value for n. We will use n=4 as an example.
# The user of this code can change the value of n to their desired value.
n_example = 4
try:
    lambda_vals = solve_for_lambda(n_example)
    print(f"For n = {n_example}, the values of λ for which the equation has no solution are:")
    for val in lambda_vals:
        print(val)
except ValueError as e:
    print(e)

# Example for n=2 to check against analytical results
# For n=2, eigenvalues of A = [[1/2, 1], [1/3, 1/2]] are 1/2 ± 1/√3.
# λ = 1 / (1/2 ± 1/√3) which are 4√3 - 6 ≈ 0.928 and -4√3 - 6 ≈ -12.928.
print("\n--- Verification for n=2 ---")
lambda_vals_n2 = solve_for_lambda(2)
print("Computed values for n=2:", lambda_vals_n2)
analytical_1 = 4 * np.sqrt(3) - 6
analytical_2 = -4 * np.sqrt(3) - 6
print("Analytical values for n=2:", np.sort([analytical_1, analytical_2]))
