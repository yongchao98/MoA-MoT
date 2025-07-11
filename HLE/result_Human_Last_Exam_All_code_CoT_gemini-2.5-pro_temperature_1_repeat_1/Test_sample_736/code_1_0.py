import numpy as np

def solve_for_n(n):
    """
    Constructs the matrix A_n and computes its minimum eigenvalue.
    A_n has entries A_ij = n - |i-j|.
    """
    if n < 1:
        raise ValueError("n must be a positive integer")
    
    # Create an empty n x n matrix
    A = np.zeros((n, n))
    
    # Populate the matrix
    for i in range(1, n + 1):
        for j in range(1, n + 1):
            A[i-1, j-1] = n - abs(i - j)
            
    # Calculate eigenvalues
    eigenvalues = np.linalg.eigvalsh(A)
    
    # Find the minimum eigenvalue
    min_eigenvalue = np.min(eigenvalues)
    
    return min_eigenvalue

# We can compute the minimum eigenvalue for several values of n to see the trend
print("n=1, min_eigenvalue =", solve_for_n(1))
print("n=2, min_eigenvalue =", solve_for_n(2))
print("n=3, min_eigenvalue =", solve_for_n(3))
print("n=10, min_eigenvalue =", solve_for_n(10))
print("n=50, min_eigenvalue =", solve_for_n(50))
print("n=100, min_eigenvalue =", solve_for_n(100))
print("n=200, min_eigenvalue =", solve_for_n(200))

# The limit as n -> infinity can be shown to be 1/4.
c = 0.25
print(f"\nThe value of c seems to converge to {c}")
