import numpy as np

def is_point_of_continuity(M):
    """
    Checks if the map M -> minimal_polynomial(M) is continuous at a given matrix M.
    This is true if and only if the matrix M is non-derogatory.
    A matrix is non-derogatory if and only if the geometric multiplicity of each of its eigenvalues is 1.

    Args:
        M (np.ndarray): A square matrix.

    Returns:
        bool: True if the map is continuous at M, False otherwise.
    """
    if M.shape[0] != M.shape[1]:
        raise ValueError("Input matrix must be square.")
    
    n = M.shape[0]
    if n == 0:
        return True
    
    try:
        eigenvalues = np.linalg.eigvals(M)
    except np.linalg.LinAlgError:
        print("Could not compute eigenvalues.")
        return False

    # Group close eigenvalues together to handle floating point inaccuracies.
    # We round them to a certain precision to find the unique ones.
    # The choice of decimal precision may be important for ill-conditioned matrices.
    unique_eigenvalues = np.unique(np.round(eigenvalues, decimals=8))

    print(f"Matrix:\n{M}")
    print(f"Eigenvalues: {eigenvalues}")
    print("-" * 30)

    for lam in unique_eigenvalues:
        # For each unique eigenvalue, compute the geometric multiplicity.
        # Geometric multiplicity = n - rank(M - lambda*I)
        I = np.identity(n)
        matrix_for_rank = M - lam * I
        rank = np.linalg.matrix_rank(matrix_for_rank)
        geom_multiplicity = n - rank
        
        print(f"For eigenvalue {lam}:")
        print(f"  Geometric multiplicity = {n} - rank(M - {lam}*I) = {n} - {rank} = {geom_multiplicity}")

        if geom_multiplicity > 1:
            print("\nResult: Found an eigenvalue with geometric multiplicity > 1.")
            print("The matrix is derogatory, so the map is discontinuous at this point.")
            return False
            
    print("\nResult: All eigenvalues have geometric multiplicity 1.")
    print("The matrix is non-derogatory, so the map is continuous at this point.")
    return True

if __name__ == '__main__':
    # Example 1: Matrix with distinct eigenvalues (non-derogatory)
    M1 = np.array([[1, 2], 
                   [3, 4]], dtype=float)
    is_point_of_continuity(M1)
    print("\n" + "="*50 + "\n")

    # Example 2: Jordan block (non-derogatory)
    M2 = np.array([[2, 1], 
                   [0, 2]], dtype=float)
    is_point_of_continuity(M2)
    print("\n" + "="*50 + "\n")

    # Example 3: Scalar matrix (derogatory)
    M3 = np.array([[2, 0], 
                   [0, 2]], dtype=float)
    is_point_of_continuity(M3)
    print("\n" + "="*50 + "\n")

    # Example 4: A 4x4 derogatory matrix
    M4 = np.array([[1, 0, 0, 0],
                   [0, 1, 0, 0],
                   [0, 0, 2, 1],
                   [0, 0, 0, 2]], dtype=float)
    is_point_of_continuity(M4)
    print("\n" + "="*50 + "\n")