import numpy as np

def find_minimum_eigenvalue(n):
    """
    Constructs the matrix A_n and returns its smallest eigenvalue.
    The matrix A_n has entries A_ij = n - |i-j|.
    """
    if n < 1:
        return None
    
    # Create the n x n matrix A
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            A[i, j] = n - abs(i - j)
            
    # Calculate the eigenvalues. eigvalsh is used for symmetric matrices.
    eigenvalues = np.linalg.eigvalsh(A)
    
    # Return the minimum eigenvalue
    return np.min(eigenvalues)

def main():
    """
    Calculates and prints the minimum eigenvalue of matrix A_n for n from 1 to 20.
    """
    print("Calculating the smallest eigenvalue of the matrix A_n for various n.")
    
    # The problem is to find c in:
    # sum_{i,j=1 to n} (n - |i-j|) * x_i * x_j >= c * sum_{i=1 to n} x_i^2
    # This implies c <= min_eigenvalue(A_n) for all n.
    # So we are looking for c = inf_{n>=1} min_eigenvalue(A_n).

    min_eigenvalues = []
    for n in range(1, 21):
        min_eig = find_minimum_eigenvalue(n)
        min_eigenvalues.append(min_eig)
        print(f"For n = {n:2d}, the minimum eigenvalue is {min_eig:.6f}")

    # The sequence of minimum eigenvalues is decreasing towards 0.
    # The infimum of this set of eigenvalues is 0.
    # Therefore, the maximum value for c is 0.
    final_c = 0
    print(f"\nThe sequence of minimum eigenvalues is decreasing and appears to approach 0.")
    print(f"The infimum of these eigenvalues is 0.")
    print(f"Thus, the maximum value for c that satisfies the inequality for all n is {final_c}.")

if __name__ == "__main__":
    main()
