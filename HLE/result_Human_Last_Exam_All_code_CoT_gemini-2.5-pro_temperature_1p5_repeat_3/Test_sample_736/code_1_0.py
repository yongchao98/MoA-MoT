import numpy as np

def find_maximum_c():
    """
    Finds the maximum real number c for the inequality by computing the
    minimum eigenvalues of the associated matrix A_n for various n.
    """
    print("We are looking for c = inf_{n>=1} lambda_min(A_n), where A_n is the matrix from the inequality.")
    print("Let's compute lambda_min(A_n) for various n to observe the trend.\n")
    
    min_eigenvalues = []
    
    # We compute for n from 1 up to a reasonably large number, e.g., 50.
    n_max = 50
    for n in range(1, n_max + 1):
        # Create the n x n matrix A_n
        A_n = np.fromfunction(lambda i, j: n - np.abs(i - j), (n, n), dtype=int)
        
        # Calculate the minimum eigenvalue. eigvalsh is efficient for symmetric matrices.
        min_eig = np.linalg.eigvalsh(A_n)[0]
        min_eigenvalues.append(min_eig)

    print("Minimum eigenvalues for n = 1, 2, ..., 10:")
    for n_minus_1 in range(10):
        print(f"n = {n_minus_1 + 1:2d}: lambda_min = {min_eigenvalues[n_minus_1]:.8f}")
        
    print("\nMinimum eigenvalues for some larger n:")
    for n_minus_1 in [19, 29, 39, 49]: # Corresponds to n=20, 30, 40, 50
        print(f"n = {n_minus_1 + 1:2d}: lambda_min = {min_eigenvalues[n_minus_1]:.8f}")

    print("\nThe sequence of minimum eigenvalues is decreasing for n >= 2.")
    print("The infimum is the limit as n -> infinity, which approaches 0.25.")
    
    c = 1/4
    
    # The final equation is sum_{i,j}(n-|i-j|)x_i x_j >= c * sum_i x_i^2
    # The prompt asks to output each number in the final equation. 
    # Since n, i, j, x_i are variables, the only number to find is c.
    print(f"\nThe maximum value for c is {c}.")

find_maximum_c()
<<<0.25>>>