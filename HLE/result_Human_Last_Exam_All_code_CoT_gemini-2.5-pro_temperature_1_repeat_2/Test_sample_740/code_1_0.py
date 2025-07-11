import numpy as np
import numpy.linalg

def find_lambda_for_no_solution(n):
    """
    Calculates the values of lambda for which the given integral equation has no solution.

    The method converts the integral equation into a matrix eigenvalue problem.
    The values of lambda for which there is no unique solution are the reciprocals
    of the eigenvalues of a matrix A, where A_jk = 1/(j + n - k).

    Args:
        n: A positive integer, the exponent in the kernel of the integral equation.

    Returns:
        A numpy array of the lambda values.
    """
    if not isinstance(n, int) or n <= 0:
        raise ValueError("n must be a positive integer.")

    # Construct the matrix A with elements A_jk = 1/(j + n - k)
    # The indices j and k run from 0 to n-1.
    j_indices, k_indices = np.ogrid[0:n, 0:n]
    A = 1.0 / (j_indices + n - k_indices)

    # Compute the eigenvalues of A.
    # The eigenvalues of A for this problem are real.
    eigenvalues_mu = numpy.linalg.eigvals(A)

    # The values of lambda for which no solution exists are the reciprocals of these eigenvalues.
    lambda_values = 1.0 / eigenvalues_mu
    
    # Sort the real parts of the values for a consistent and clean output.
    return np.sort(np.real(lambda_values))

# The user can specify the integer 'n'.
# We choose n=2 as an example, for which the analytical solution is known.
n = 2
lambda_vals = find_lambda_for_no_solution(n)

print(f"For n = {n}, the values of λ for which the equation has no solution are:")
for val in lambda_vals:
    print(f"{val:.8f}")

# For n=4, we can find the values by changing the value of n.
n = 4
lambda_vals_n4 = find_lambda_for_no_solution(n)
print(f"\nFor n = {n}, the values of λ for which the equation has no solution are:")
for val in lambda_vals_n4:
    print(f"{val:.8f}")
