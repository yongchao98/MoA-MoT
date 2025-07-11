import numpy as np

def check_continuity_point(M):
    """
    Checks if a matrix M is a point of continuity for the minimal polynomial map.
    This is equivalent to checking if M is non-derogatory.

    A matrix is non-derogatory if and only if the geometric multiplicity of
    each of its eigenvalues is 1.
    Geometric multiplicity of eigenvalue lambda = n - rank(M - lambda*I).
    """
    print(f"--- Checking matrix ---\n{M}\n")
    try:
        if M.shape[0] != M.shape[1]:
            print("Matrix must be square.")
            return False
        n = M.shape[0]
        
        # Calculate eigenvalues. Use a tolerance for uniqueness check.
        eigenvalues = np.linalg.eigvals(M)
        
        # Group close eigenvalues to find unique ones numerically
        tol = 1e-8
        eigenvalues.sort()
        unique_eigenvalues = []
        if len(eigenvalues) > 0:
            unique_eigenvalues.append(eigenvalues[0])
            for i in range(1, len(eigenvalues)):
                if abs(eigenvalues[i] - unique_eigenvalues[-1]) > tol:
                    unique_eigenvalues.append(eigenvalues[i])
        
        print(f"Unique eigenvalues found: {[np.round(e, 5) for e in unique_eigenvalues]}")

        is_non_derogatory = True
        for lam in unique_eigenvalues:
            # Form the matrix M - lambda*I
            mat = M - lam * np.identity(n)
            
            # Calculate rank
            rank = np.linalg.matrix_rank(mat, tol=tol)
            
            # Calculate geometric multiplicity
            geom_mult = n - rank
            
            # Output the equation for geometric multiplicity
            print(f"For eigenvalue lambda ~ {np.round(lam.real, 5) if np.abs(lam.imag) < tol else np.round(lam, 5)}:")
            print(f"  Geometric multiplicity = n - rank(M - lambda*I) = {n} - {rank} = {geom_mult}")

            if geom_mult > 1:
                is_non_derogatory = False
                print(f"  Multiplicity > 1. Matrix is derogatory.")
                break
            else:
                print("  Multiplicity is 1. OK.")

        print("-" * 23)
        if is_non_derogatory:
            print("Result: All eigenvalues have geometric multiplicity 1.")
            print("The matrix is NON-DEROGATORY and thus a POINT OF CONTINUITY.")
        else:
            print("Result: At least one eigenvalue has geometric multiplicity > 1.")
            print("The matrix is DEROGATORY and thus a POINT OF DISCONTINUITY.")
        print("\n")
        return is_non_derogatory

    except Exception as e:
        print(f"An error occurred: {e}")
        return False

if __name__ == '__main__':
    # Example 1: A derogatory matrix (a point of discontinuity)
    # This is a scalar matrix. Eigenvalue is 2, with algebraic multiplicity 2.
    # The eigenspace is the whole R^2, so geometric multiplicity is 2.
    M1 = np.array([[2, 0], 
                   [0, 2]], dtype=float)
    check_continuity_point(M1)

    # Example 2: A non-derogatory matrix (a point of continuity)
    # This matrix has a single eigenvalue 2, with algebraic multiplicity 2.
    # The eigenspace is spanned by (1, 0), so geometric multiplicity is 1.
    M2 = np.array([[2, 1], 
                   [0, 2]], dtype=float)
    check_continuity_point(M2)

    # Example 3: A non-derogatory matrix with distinct eigenvalues (a point of continuity)
    M3 = np.array([[1, 2, 3], 
                   [4, 5, 6], 
                   [7, 8, 10]], dtype=float)
    check_continuity_point(M3)

    # Example 4: A 3x3 derogatory matrix (a point of discontinuity)
    # This matrix has eigenvalues 2, 2, 3.
    # For lambda=2, M-2I has rank 1, so geometric multiplicity is 3-1=2.
    M4 = np.array([[2, 0, 0],
                   [0, 2, 1],
                   [0, 0, 3]], dtype=float)
    check_continuity_point(M4)
