import numpy as np

def find_min_eigenvalue(n):
    """
    Constructs the n x n matrix A and calculates its minimum eigenvalue.
    A_ij = n - |i - j|
    """
    if n <= 0:
        return float('inf')
    
    # Create the matrix A_n
    A = np.fromfunction(lambda i, j: n - np.abs(i - j), (n, n), dtype=np.float64)
    
    # Calculate eigenvalues. For a symmetric matrix, eigvalsh is efficient.
    eigenvalues = np.linalg.eigvalsh(A)
    
    # Return the minimum eigenvalue
    return np.min(eigenvalues)

def main():
    """
    Calculates the minimum eigenvalue for n=1 to 20 to show the trend,
    and then for a large n to approximate c.
    """
    print("Minimum eigenvalues for n=1 to 20:")
    for n in range(1, 21):
        min_eig = find_min_eigenvalue(n)
        print(f"For n = {n:2d}, the minimum eigenvalue is: {min_eig:.15f}")
        
    print("\nApproximating c with a larger n:")
    n_large = 200
    c_approx = find_min_eigenvalue(n_large)
    print(f"For n = {n_large}, the minimum eigenvalue is: {c_approx:.15f}")
    print("\nThe maximum value for c is the limit of this decreasing sequence.")
    print(f"An approximation for c is {c_approx:.15f}")

if __name__ == "__main__":
    main()
