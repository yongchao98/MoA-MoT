import numpy as np
import cmath

def find_lambda_for_no_solution(n):
    """
    Calculates the values of lambda for which the given integral equation has no solution.
    
    The problem is reduced to finding the reciprocals of the eigenvalues of a matrix A,
    where A_jk = 1 / (j + n - k) for j,k in {0, 1, ..., n-1}.
    
    Args:
        n (int): The parameter n from the integral equation.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.")
        return

    # Construct the n x n matrix A
    A = np.zeros((n, n), dtype=float)
    for j in range(n):
        for k in range(n):
            A[j, k] = 1.0 / (j + n - k)
            
    # Calculate the eigenvalues of A. They can be real or complex.
    try:
        eigenvalues = np.linalg.eigvals(A)
    except np.linalg.LinAlgError:
        print("Error: Eigenvalue computation failed.")
        return
        
    # The values of lambda are the reciprocals of the eigenvalues of A.
    # We filter out any potential zero eigenvalues to avoid division by zero,
    # although for this specific matrix A, the eigenvalues are non-zero.
    lambda_values = [1.0 / mu for mu in eigenvalues if abs(mu) > 1e-12]

    print(f"For n = {n}, the values of lambda for which the equation has no solution are:")
    for val in lambda_values:
        # We check if the imaginary part is negligible and print accordingly.
        if abs(val.imag) < 1e-9:
            print(f"{val.real:.8f}")
        else:
            print(f"{val.real:.8f} + {val.imag:.8f}j")

# We will solve for the case n=2 as a demonstration.
# The analytical solution for n=2 is λ = -6 ± 4√3.
# Let's compare the numerical results.
# 4√3 is approximately 4 * 1.73205 = 6.9282
# λ ≈ -6 + 6.9282 = 0.9282
# λ ≈ -6 - 6.9282 = -12.9282
n_example = 2
find_lambda_for_no_solution(n_example)