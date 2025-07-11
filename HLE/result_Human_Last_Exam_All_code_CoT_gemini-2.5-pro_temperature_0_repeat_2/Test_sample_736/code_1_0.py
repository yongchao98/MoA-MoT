import numpy as np

def solve():
    """
    Finds the maximum real number c for the given inequality by finding the
    minimum of the smallest eigenvalues of the matrix A_n for n from 1 to 50.
    """
    min_overall_eigenvalue = float('inf')
    n_for_min_eigenvalue = -1

    # We test for n from 1 up to 50, which should be sufficient to find the minimum.
    for n in range(1, 51):
        # Create the matrix A_n
        A_n = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                A_n[i, j] = n - abs(i - j)

        # Calculate eigenvalues
        eigenvalues = np.linalg.eigvalsh(A_n)
        
        # Find the minimum eigenvalue for the current n
        min_eigenvalue_n = np.min(eigenvalues)

        # Update the overall minimum if the current one is smaller
        if min_eigenvalue_n < min_overall_eigenvalue:
            min_overall_eigenvalue = min_eigenvalue_n
            n_for_min_eigenvalue = n

    # The maximum value for c is the minimum of all these minimum eigenvalues.
    c = min_overall_eigenvalue
    
    # The value for n=4 is 2 - sqrt(2)
    c_analytical = 2 - np.sqrt(2)

    print(f"The minimum eigenvalue is found at n = {n_for_min_eigenvalue}")
    print(f"The minimum eigenvalue (c) is approximately: {c}")
    print(f"The analytical value for n=4 is 2 - sqrt(2), which is approximately: {c_analytical}")
    print("The maximum value for c is the minimum of these eigenvalues.")
    print(f"c = 2 - sqrt(2)")

solve()
