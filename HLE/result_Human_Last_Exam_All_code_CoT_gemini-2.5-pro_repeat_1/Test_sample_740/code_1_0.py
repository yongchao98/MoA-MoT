import numpy as np

def solve_for_lambda(n):
    """
    Computes the values of lambda for which the integral equation has no solution.

    Args:
        n (int): The exponent in the integral equation kernel.

    Returns:
        numpy.ndarray: An array of lambda values.
    """
    if n <= 0:
        raise ValueError("n must be a positive integer")

    # Construct the matrix A
    # A_jk = 1 / (j + n - k) where j is row index, k is column index
    A = np.zeros((n, n))
    for j in range(n):
        for k in range(n):
            A[j, k] = 1.0 / (j + n - k)

    # Find the eigenvalues of A
    eigenvalues_mu = np.linalg.eigvals(A)

    # The values of lambda are the reciprocals of the eigenvalues of A
    # We filter out potential division by zero, though eigenvalues are non-zero for this matrix
    lambda_values = np.array([1/mu for mu in eigenvalues_mu if mu != 0])
    
    return lambda_values

# Let's solve for a specific case, e.g., n = 3.
n_val = 3
try:
    lambdas = solve_for_lambda(n_val)
    print(f"For n = {n_val}, the values of lambda for which the equation has no solution are:")
    for i, lam in enumerate(lambdas):
        # "Remember in the final code you still need to output each number in the final equation!"
        # Printing the values in the format "lambda_i = value"
        print(f"λ_{i+1} = {lam:.8f}")
except ValueError as e:
    print(e)
