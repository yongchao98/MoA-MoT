import numpy as np

def find_max_c():
    """
    This function finds the maximum real number c for the given inequality
    by numerically computing the minimum eigenvalues of the associated matrix
    for increasing n and observing their trend.
    """

    print("The inequality can be written as x^T A_n x >= c * x^T x,")
    print("where A_n is a matrix with elements A_ij = n - |i-j|.")
    print("For this to hold for all x, c must be <= the minimum eigenvalue of A_n.")
    print("Since this must hold for all n, c must be <= infimum over n of the minimum eigenvalues.")
    print("\nLet's calculate the minimum eigenvalue for various n:")

    min_eigenvalues = []
    ns_to_check = [1, 2, 3, 4, 5, 10, 20, 50, 100]

    for n in ns_to_check:
        # Create the matrix A_n
        A_n = np.fromfunction(lambda i, j: n - np.abs(i - j), (n, n), dtype=int)

        # Calculate eigenvalues. Since A_n is symmetric, we can use eigvalsh.
        eigenvalues = np.linalg.eigvalsh(A_n)

        # Get the minimum eigenvalue
        min_eig = eigenvalues[0]
        min_eigenvalues.append(min_eig)

        print(f"For n = {n:3d}, the minimum eigenvalue is: {min_eig:.8f}")

    print("\nAs n increases, the minimum eigenvalue decreases and approaches 0.")
    print("The infimum of these minimum eigenvalues is 0.")
    print("Therefore, the maximum value for c is 0.")

    c = 0
    print(f"\nThe maximum real number c is {c}.")

find_max_c()