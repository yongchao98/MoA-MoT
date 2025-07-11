import numpy as np

def find_lambda_values_with_no_solution(n):
    """
    Calculates the values of λ for which the given integral equation has no solution.

    The problem is reduced to finding the reciprocals of the eigenvalues of a
    specific matrix M.

    Args:
        n (int): The integer parameter n from the integral equation.

    Returns:
        numpy.ndarray: An array containing the complex values of λ.
    """
    # Based on the analysis, the critical values of λ are the reciprocals
    # of the eigenvalues of the matrix M, where M_ij = 1/(n + j - i).
    # i denotes the row index (from 0 to n-1) and j the column index (from 0 to n-1).
    M = np.array([[1.0 / (n + j - i) for j in range(n)] for i in range(n)])

    # Calculate the eigenvalues of the matrix M.
    try:
        eigenvalues = np.linalg.eigvals(M)
    except np.linalg.LinAlgError:
        print(f"Error: Eigenvalue computation failed for n={n}. The matrix might be singular or ill-conditioned.")
        return np.array([])

    # The values of λ for which no solution exists are the reciprocals of these eigenvalues.
    # We must handle the case of zero eigenvalues, which would correspond to λ -> infinity.
    # We filter out numerically small eigenvalues to avoid division by zero.
    # Based on symbolic checks for small n, M is non-singular.
    stable_eigenvalues = eigenvalues[np.abs(eigenvalues) > 1e-12]

    if len(stable_eigenvalues) < n:
        print(f"Warning: For n={n}, the matrix M appears to be singular (has zero or near-zero eigenvalues).")

    # λ = 1 / μ
    lambdas = 1.0 / stable_eigenvalues
    
    # Sort the values for consistent output, by real part then imaginary part.
    lambdas_sorted = sorted(lambdas, key=lambda x: (x.real, x.imag))

    return np.array(lambdas_sorted)

def main():
    """
    Demonstrates the calculation for n = 1, 2, and 3 and prints the results.
    """
    print("This script calculates the values of λ for which the equation has no solutions.")
    
    for n_val in [1, 2, 3]:
        print("-" * 30)
        print(f"For the case n = {n_val}:")
        
        lambdas = find_lambda_values_with_no_solution(n_val)
        
        if lambdas.size == 0:
            print("No finite values of λ found.")
            continue
            
        print("The values of λ are:")
        for lam in lambdas:
            # Print the values in a readable format.
            if abs(lam.imag) < 1e-9:
                # Treat as a real number
                print(f"λ = {lam.real:.8f}")
            else:
                # Print as a complex number
                print(f"λ = {lam.real:.8f} + {lam.imag:.8f}j")

main()