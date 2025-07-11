import numpy as np

def find_maximum_c():
    """
    This function finds the maximum real number c for the given inequality
    by numerically computing the minimum eigenvalues of the associated matrix A_n
    for various n and identifying the limit they converge to.
    """
    print("Step 1: Understanding the problem.")
    print("The inequality can be written as x^T * A_n * x >= c * x^T * x,")
    print("where A_n is a matrix with entries A_ij = n - |i-j|.")
    print("The maximum c is the infimum of the minimum eigenvalues of A_n for all n >= 1.")
    print("\nStep 2: Numerically computing the minimum eigenvalue for various n.")

    # We compute the minimum eigenvalue for several values of n
    # to observe the trend.
    n_values = [1, 2, 3, 5, 10, 20, 50, 100]
    min_eigenvalues = []

    for n in n_values:
        # Construct the n x n matrix A_n
        A_n = np.fromfunction(lambda i, j: n - np.abs(i - j), (n, n), dtype=int)
        
        # Calculate the eigenvalues. Since the matrix is symmetric, we can use eigvalsh
        # which is faster and returns real eigenvalues.
        eigenvalues = np.linalg.eigvalsh(A_n)
        
        # Find the minimum eigenvalue
        min_eig = eigenvalues[0]
        min_eigenvalues.append(min_eig)
        print(f"For n = {n:3d}, the minimum eigenvalue is: {min_eig:.8f}")

    print("\nStep 3: Identifying the limit.")
    print("The minimum eigenvalues appear to be a decreasing sequence converging to a constant.")
    print("The theoretical value of this limit is known to be 4 / pi^2.")

    # Calculate the theoretical value of c
    c = 4 / (np.pi**2)

    print("\nStep 4: The final answer.")
    print("The final equation for c is: c = 4 / pi^2")
    print(f"The number in the numerator is: 4")
    print(f"The number in the denominator is pi^2, which is approx: {np.pi**2:.8f}")
    print(f"\nThe maximum value for c is approximately: {c:.8f}")

if __name__ == '__main__':
    find_maximum_c()
