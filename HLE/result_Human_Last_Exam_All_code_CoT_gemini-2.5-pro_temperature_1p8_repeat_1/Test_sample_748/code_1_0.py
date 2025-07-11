import numpy as np
from numpy.linalg import matrix_power, matrix_rank

def get_minimal_poly_degree(M: np.ndarray) -> int:
    """
    Calculates the degree of the minimal polynomial for a square matrix M.

    Args:
        M: A numpy array representing the square matrix.

    Returns:
        The degree of the minimal polynomial of M.
    """
    n = M.shape[0]
    if n == 0:
        return 0

    # Tolerance for floating-point comparisons
    tol = 1e-8

    # Step 1: Compute eigenvalues
    eigenvalues = np.linalg.eigvals(M)
    
    # Step 2: Find distinct eigenvalues by grouping close ones
    distinct_eigenvalues = []
    # Create a copy to sort and modify
    temp_eigs = sorted(eigenvalues.tolist(), key=lambda x: (np.real(x), np.imag(x)))
    
    while temp_eigs:
        mu = temp_eigs.pop(0)
        distinct_eigenvalues.append(mu)
        # Remove other eigenvalues that are numerically close to mu
        temp_eigs = [eig for eig in temp_eigs if np.abs(eig - mu) > tol]

    min_poly_degree = 0
    
    # Step 3: For each distinct eigenvalue, find its index m_j
    for mu in distinct_eigenvalues:
        A = M - mu * np.identity(n)
        
        # Find the smallest integer m >= 1 s.t. rank(A^m) == rank(A^(m+1))
        # This m is the index of the eigenvalue mu.
        m = 1
        while m <= n:
            Am = matrix_power(A, m)
            rank_m = matrix_rank(Am, tol=tol)
            
            A_m_plus_1 = matrix_power(A, m + 1)
            rank_m_plus_1 = matrix_rank(A_m_plus_1, tol=tol)
            
            if rank_m == rank_m_plus_1:
                min_poly_degree += m
                break
            m += 1
        else:
            # This case should ideally not be reached in exact arithmetic,
            # as m <= algebraic multiplicity <= n.
            # If it is, it might indicate numerical instability.
            # Let's assume the index is n, which ensures we classify it as non-derogatory.
            min_poly_degree += n
            break # No need to check other eigenvalues if one has index n.
            
    return min_poly_degree

def check_continuity_point(M: np.ndarray):
    """
    Checks if a matrix M is a point of continuity for the minimal polynomial map.
    """
    n = M.shape[0]
    print("="*30)
    print(f"Testing Matrix:\n{M}")
    print(f"Size of matrix, n = {n}")

    if n == 0:
        print("The 0x0 matrix is trivially a point of continuity.")
        return

    deg = get_minimal_poly_degree(M)
    print(f"Degree of minimal polynomial = {deg}")
    
    is_continuous = (deg == n)
    
    print(f"Is M a point of continuity? -> {is_continuous}")
    if is_continuous:
        print("(The matrix is non-derogatory)")
    else:
        print("(The matrix is derogatory)")
    print("="*30 + "\n")


# --- Example Matrices ---

# Example 1: A non-derogatory matrix (companion matrix)
# Minimal polynomial = Characteristic polynomial = x^3 - x^2 - x - 1
M1 = np.array([[0, 0, 1],
               [1, 0, 1],
               [0, 1, 1]], dtype=float)
check_continuity_point(M1)

# Example 2: A derogatory matrix (scalar matrix, 2*I)
# Minimal polynomial is x-2, degree 1. Characteristic is (x-2)^3, degree 3.
M2 = np.array([[2, 0, 0],
               [0, 2, 0],
               [0, 0, 2]], dtype=float)
check_continuity_point(M2)

# Example 3: A derogatory matrix with a 2x2 Jordan block.
# Minimal poly is (x-3)^2, degree 2. Characteristic is (x-3)^3, degree 3.
M3 = np.array([[3, 1, 0],
               [0, 3, 0],
               [0, 0, 3]], dtype=float)
check_continuity_point(M3)

# Example 4: A 2x2 matrix with distinct eigenvalues, which is always non-derogatory.
M4 = np.array([[1, 2],
               [3, 4]], dtype=float)
check_continuity_point(M4)
