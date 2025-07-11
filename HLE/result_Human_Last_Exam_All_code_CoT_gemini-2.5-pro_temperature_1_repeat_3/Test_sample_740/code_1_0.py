import numpy as np

def solve_for_lambda(n):
    """
    Calculates the values of λ for which the given integral equation has no solutions.

    Args:
        n (int): The exponent in the kernel of the integral equation. 
                 It must be a positive integer.

    The method involves finding the eigenvalues of a specific n x n matrix M,
    where M_mj = 1 / (n + j - m). The values of λ for which no solution exists
    are the reciprocals of these eigenvalues.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.")
        return

    # Construct the n x n matrix M
    M = np.array([[1.0 / (n + j - m) for j in range(n)] for m in range(n)])

    # Find the eigenvalues of M
    try:
        eigenvalues = np.linalg.eigvals(M)
    except np.linalg.LinAlgError as e:
        print(f"Error computing eigenvalues for n={n}: {e}")
        return

    # The values of λ are the reciprocals of the eigenvalues of M.
    # We sort the results for consistent ordering.
    lambda_values = sorted([1.0 / mu for mu in eigenvalues if abs(mu) > 1e-15])

    print(f"For n = {n}, the values of λ for which the equation has no solution are:")
    for val in lambda_values:
        print(val)
    
    return lambda_values

if __name__ == '__main__':
    # The problem is defined for a general positive integer 'n'.
    # You can change the value of n_param below to solve for a different case.
    # For this example, we will solve for n = 2.
    n_param = 2
    
    # Run the solver
    solve_for_lambda(n_param)