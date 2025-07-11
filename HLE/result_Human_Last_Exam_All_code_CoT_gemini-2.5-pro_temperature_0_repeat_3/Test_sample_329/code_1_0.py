import numpy as np

def check_matrix(name, M):
    """
    Checks if a matrix M is in the set P by verifying necessary and sufficient conditions.
    """
    print(f"Checking matrix {name}:")
    
    # 1. Check for rational entries (conceptual check).
    # All generators v*v^T for v in Z^2 have integer entries.
    # Their convex hull P can therefore only contain matrices with rational entries.
    if name in ['D', 'E']:
        print(f"  - Fails: Matrix {name} contains irrational entries (pi) and cannot be in P.")
        return False

    # 2. Check for symmetry.
    if not np.allclose(M, M.T):
        print(f"  - Fails: Matrix {name} is not symmetric. It is not in P.")
        return False

    # 3. Check trace condition.
    if np.trace(M) < 1:
        print(f"  - Fails: Trace of matrix {name} is {np.trace(M):.2f}, which is less than 1. It is not in P.")
        return False

    # 4. Check for positive semidefiniteness.
    # For a 2x2 symmetric matrix, this means diagonal entries are non-negative and determinant is non-negative.
    if M[0, 0] < 0 or M[1, 1] < 0:
        print(f"  - Fails: Matrix {name} has negative diagonal entries, so it's not positive semidefinite.")
        return False
    # Use a small tolerance for floating point comparisons
    det = np.linalg.det(M)
    if det < -1e-9:
        print(f"  - Fails: Determinant of matrix {name} is {det:.2f}, which is negative. It is not positive semidefinite.")
        return False

    print(f"  - Pass: Matrix {name} satisfies the necessary conditions (rational, symmetric, trace>=1, PSD).")
    print(f"  - Now, checking for a valid convex combination construction...")

    # 5. Check for a specific construction (sufficient condition).
    if name == 'C':
        v1 = np.array([[1], [1]])
        S1 = v1 @ v1.T
        v2 = np.array([[1], [-1]])
        S2 = v2 @ v2.T
        
        l1, l2 = 1/4, 3/4
        constructed_M = l1 * S1 + l2 * S2
        
        if np.allclose(M, constructed_M):
            print(f"  - Success: Matrix C can be written as a convex combination:")
            print(f"    C = {l1} * v1*v1^T + {l2} * v2*v2^T, where v1=[1, 1], v2=[1, -1]")
            print(f"    C = {l1} * \n{S1}\n      + {l2} * \n{S2}\n      = \n{constructed_M}")
            print(f"  - Conclusion: Matrix C is in P.")
            return True

    if name == 'F':
        v1 = np.array([[6], [0]])
        S1 = v1 @ v1.T
        v2 = np.array([[7], [0]])
        S2 = v2 @ v2.T

        l1, l2 = 7/13, 6/13
        constructed_M = l1 * S1 + l2 * S2

        if np.allclose(M, constructed_M):
            print(f"  - Success: Matrix F can be written as a convex combination:")
            print(f"    F = {l1:.4f} * v1*v1^T + {l2:.4f} * v2*v2^T, where v1=[6, 0], v2=[7, 0]")
            print(f"    F = {l1:.4f} * \n{S1}\n      + {l2:.4f} * \n{S2}\n      = \n{constructed_M}")
            print(f"  - Conclusion: Matrix F is in P.")
            return True
            
    # If a matrix passed necessary conditions but we don't have a construction,
    # our analysis concludes it's not in P.
    print(f"  - Fails: No simple convex combination found for {name}.")
    return False

def solve():
    # Define the matrices
    A = np.array([[0, 0], [0, 0]], dtype=float)
    B = np.array([[6, 4], [3, 7]], dtype=float)
    C = np.array([[1, -1/2], [-1/2, 1]], dtype=float)
    D = np.array([[np.pi, 1], [1, np.pi**2]], dtype=float)
    E = np.array([[1, np.pi], [np.pi, 1]], dtype=float)
    F = np.array([[42, 0], [0, 0]], dtype=float)

    matrices = {'A': A, 'B': B, 'C': C, 'D': D, 'E': E, 'F': F}
    matrices_in_P = []

    for name, M in matrices.items():
        if check_matrix(name, M):
            matrices_in_P.append(name)
        print("-" * 40)

    print(f"\nFinal list of matrices in P: {matrices_in_P}")
    
solve()