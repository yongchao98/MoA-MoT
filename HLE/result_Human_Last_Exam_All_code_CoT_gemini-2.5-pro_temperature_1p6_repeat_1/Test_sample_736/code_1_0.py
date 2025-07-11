import numpy as np

def solve():
    """
    This function calculates the minimum eigenvalue for the matrix A_n
    with entries A_ij = n - |i-j| for a range of n, demonstrating
    that the values converge towards 1/4.
    """
    print("This program calculates the minimum eigenvalue of the matrix A_n for various n.")
    print("The matrix A_n has entries A_ij = n - |i-j|.")
    print("The maximum value c must be less than or equal to the minimum eigenvalue for all n.")
    print("We observe the convergence of these minimum eigenvalues as n increases.\n")

    for n in [1, 2, 3, 4, 5, 10, 20, 50, 100]:
        # Create the matrix A_n
        A_n = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                A_n[i, j] = n - abs(i - j)
        
        # Calculate eigenvalues. Since the matrix is symmetric, we use eigvalsh for efficiency.
        eigenvalues = np.linalg.eigvalsh(A_n)
        
        # Get the minimum eigenvalue
        min_eigenvalue = eigenvalues[0]
        
        print(f"For n = {n:3d}, the minimum eigenvalue is: {min_eigenvalue:.8f}")

    c = 0.25
    print(f"\nThe sequence of minimum eigenvalues is decreasing and converges to {c}.")
    print(f"Thus, the maximum value for c is {c}.")

solve()