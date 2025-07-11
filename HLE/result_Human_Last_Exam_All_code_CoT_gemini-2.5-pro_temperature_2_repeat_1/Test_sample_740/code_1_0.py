import numpy as np
import sys

def solve_for_lambda(n):
    """
    Calculates the values of lambda for which the integral equation has no solution.

    Args:
        n (int): The exponent 'n' in the integral equation. Must be a positive integer.

    Returns:
        numpy.ndarray: An array of lambda values.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: 'n' must be a positive integer.")
        return None

    # Construct the matrix M based on the derivation.
    # M_ij = 1 / (n + i - j)
    M = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            M[i, j] = 1.0 / (n + i - j)

    # The critical values of lambda are the reciprocals of the eigenvalues of M.
    try:
        eigenvalues = np.linalg.eigvals(M)
    except np.linalg.LinAlgError:
        print(f"Error: Eigenvalue computation failed for n={n}.")
        return None
        
    # Filter out eigenvalues that are very close to zero to avoid division by zero.
    non_zero_eigenvalues = eigenvalues[np.abs(eigenvalues) > 1e-9]
    
    if len(non_zero_eigenvalues) < len(eigenvalues):
        print("Warning: Matrix M may be singular (has zero eigenvalues).")
    
    lambdas = 1.0 / non_zero_eigenvalues
    
    # Sort the real part for consistent output order
    lambdas = lambdas[np.argsort(np.real(lambdas))]
    
    return lambdas

def main():
    """
    Main function to run the calculation and print the results.
    """
    # You can change the value of n here.
    try:
        n_param = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    except ValueError:
        print("Invalid input. Please provide an integer value for n.")
        n_param = 4

    print(f"Calculating values of lambda for n = {n_param}\n")
    
    lambda_values = solve_for_lambda(n_param)

    if lambda_values is not None:
        print("The values of λ for which the equation has no solutions are:")
        # The instruction "output each number in the final equation" is interpreted
        # as clearly printing each value of lambda in its own equation form.
        for val in lambda_values:
            # np.linalg.eigvals can return complex numbers for non-symmetric matrices,
            # even if eigenvalues are theoretically real. We print them accordingly.
            if np.imag(val) == 0:
                print(f"λ = {np.real(val):.8f}")
            else:
                print(f"λ = {val:.8f}")

if __name__ == '__main__':
    main()
