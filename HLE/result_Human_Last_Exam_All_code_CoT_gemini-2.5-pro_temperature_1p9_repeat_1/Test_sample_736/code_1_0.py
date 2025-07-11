import numpy as np
import math

def find_c(n_large=1000):
    """
    This function finds the maximum value of c by calculating the minimum eigenvalue
    of the matrix A_n for a large n, which approximates the limit as n -> infinity.
    """
    # Create the n x n matrix A_n where A_ij = n - |i-j|
    # Note: for i,j starting at 0 for array indexing, this is n - |(i+1)-(j+1)| = n - |i-j|
    A_n = np.fromfunction(lambda i, j: n_large - np.abs(i - j), (n_large, n_large))
    
    # Calculate the eigenvalues. Since A_n is symmetric, we can use eigvalsh for efficiency.
    eigenvalues = np.linalg.eigvalsh(A_n)
    
    # The minimum eigenvalue approximates c
    c_approx = np.min(eigenvalues)
    
    pi = math.pi
    c_exact_val = 2 / (pi**2)

    print(f"The numerical approximation for c with n={n_large} is: {c_approx:.8f}")
    print(f"The theoretical value for c is 2 / pi^2.")
    print(f"Using pi ≈ {pi:.8f}, we have the equation:")
    print(f"c = 2 / ({pi:.8f}^2) = 2 / {pi**2:.8f} = {c_exact_val:.8f}")

if __name__ == '__main__':
    find_c()
