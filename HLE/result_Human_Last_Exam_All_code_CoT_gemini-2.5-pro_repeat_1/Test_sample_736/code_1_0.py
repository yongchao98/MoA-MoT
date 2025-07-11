import numpy as np

def solve_for_c(n):
    """
    Calculates the minimum eigenvalue for the matrix A_n.
    A_ij = n - |i-j|
    """
    # Create an empty n x n matrix
    A = np.zeros((n, n))
    
    # Populate the matrix based on the formula
    for i in range(n):
        for j in range(n):
            A[i, j] = n - abs(i - j)
            
    # Calculate the eigenvalues. Since A is symmetric, we use eigvalsh for efficiency.
    eigenvalues = np.linalg.eigvalsh(A)
    
    # The minimum eigenvalue is the first one in the sorted list.
    return eigenvalues[0]

# For small n, the eigenvalues can be found analytically or easily computed.
# n=1 -> 1.0
# n=2 -> 1.0
# n=3 -> (7 - sqrt(33))/2 ~ 0.627
# n=4 -> ~0.420
# The sequence seems to decrease, so the infimum will be the limit.
# We approximate this limit by taking a large value for n.
n_large = 4000
c = solve_for_c(n_large)

# The problem is to find c in the inequality. Let's output the found value.
print(f"The maximum real number c is approximately: {c}")
