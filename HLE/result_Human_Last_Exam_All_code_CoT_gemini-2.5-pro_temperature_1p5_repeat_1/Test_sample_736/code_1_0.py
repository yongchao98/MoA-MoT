import numpy as np

def find_min_eigenvalue(n):
    """
    Constructs the matrix A_n and computes its minimum eigenvalue.
    A_ij = n - |i-j|
    """
    if n == 0:
        return None
    matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            matrix[i, j] = n - abs(i - j)
            
    # For symmetric/Hermitian matrices, eigvalsh is preferred
    eigenvalues = np.linalg.eigvalsh(matrix)
    return np.min(eigenvalues)

def main():
    print("n | min_eigenvalue(A_n)")
    print("------------------------")
    for n in range(1, 21):
        min_eig = find_min_eigenvalue(n)
        print(f"{n:2d} | {min_eig:.8f}")
        
    print("\nObserving the trend for larger n:")
    for n in [50, 100, 200]:
        min_eig = find_min_eigenvalue(n)
        print(f"{n:3d} | {min_eig:.8f}")

if __name__ == '__main__':
    main()
