import numpy as np

def find_lambda_for_no_solution(n):
    """
    Calculates the values of lambda for which the integral equation has no solution.

    Args:
        n (int): The exponent in the integral equation kernel. Must be a positive integer.

    Returns:
        A sorted list of floats representing the values of lambda.
    """
    if not isinstance(n, int) or n < 1:
        raise ValueError("n must be a positive integer.")

    # Create the n x n matrix H
    # H_pj = 1 / (n + j - p)
    # p is the row index, j is the column index
    H = np.zeros((n, n))
    for p in range(n):
        for j in range(n):
            H[p, j] = 1.0 / (n + j - p)

    # The values of 1/lambda are the eigenvalues of H.
    # So, lambda = 1 / eigenvalue.
    # The matrix H is the product of an anti-diagonal matrix and the Hilbert
    # matrix, ensuring its eigenvalues are real and non-zero.
    try:
        eigenvalues = np.linalg.eigvals(H)
    except np.linalg.LinAlgError:
        print(f"Eigenvalue computation failed for n={n}.")
        return []

    # Calculate lambda values by taking the reciprocal of the eigenvalues
    lambda_values = 1.0 / eigenvalues

    # Sort the values for consistent output
    lambda_values.sort()

    return lambda_values

if __name__ == "__main__":
    # We solve for the case n=2, as discussed in the plan.
    # The user can change this value to solve for other cases.
    n_val = 2
    
    # For n=2, the exact symbolic solution is lambda = -6 ± 4*sqrt(3)
    # Symbolic values:
    val1_symbolic = -6 + 4 * np.sqrt(3)  # approximately 0.9282
    val2_symbolic = -6 - 4 * np.sqrt(3)  # approximately -12.9282
    
    # Calculate numerically using our function
    lambda_vals_numeric = find_lambda_for_no_solution(n_val)

    print(f"For n = {n_val}, the values of lambda for which the equation has no solution are:")
    # We output the numeric values which should match the symbolic ones.
    print(np.round(lambda_vals_numeric, 8))