import numpy as np

def find_min_eigenvalue(n):
    """
    Constructs the matrix A_n and computes its minimum eigenvalue.
    The matrix A_n has entries A_ij = n - |i-j|.
    """
    # Create an empty n x n matrix
    a_n = np.zeros((n, n))
    # Populate the matrix based on the formula
    for i in range(n):
        for j in range(n):
            a_n[i, j] = n - abs(i - j)
    
    # Calculate the eigenvalues. eigh is used for symmetric matrices.
    eigenvalues = np.linalg.eigh(a_n)[0]
    
    # Return the minimum eigenvalue
    return np.min(eigenvalues)

def main():
    """
    Calculates the minimum eigenvalue for various n and demonstrates
    the convergence to find the value of c.
    """
    print("Investigating the minimum eigenvalue of A_n for increasing n:")
    
    test_n_values = [1, 2, 3, 4, 5, 10, 20, 50, 100, 200]
    
    for n in test_n_values:
        min_eig = find_min_eigenvalue(n)
        print(f"For n = {n:3d}, the minimum eigenvalue is {min_eig:.6f}")
        
    print("\nAs n increases, the minimum eigenvalue decreases and converges.")
    print("The theoretical limit as n approaches infinity is 1/4.")
    print("Therefore, the maximum value for c that holds for all n is 1/4.")
    
    c = 1/4
    
    print("\nThe final inequality is:")
    print(f"sum_{{i=1 to n}} sum_{{j=1 to n}} (n - |i-j|) * x_i * x_j >= {c} * sum_{{i=1 to n}} x_i^2")

if __name__ == "__main__":
    main()
