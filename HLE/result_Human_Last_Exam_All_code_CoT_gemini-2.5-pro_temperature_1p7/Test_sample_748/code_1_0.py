import numpy as np

def is_point_of_continuity(M):
    """
    Checks if a matrix M is a point of continuity for the minimal polynomial map.
    This is equivalent to checking if M is a non-derogatory matrix.

    A matrix M is non-derogatory if and only if for each of its distinct eigenvalues lambda,
    the geometric multiplicity is 1. This is equivalent to rank(M - lambda*I) = n-1.

    Args:
        M (np.ndarray): An n x n matrix.

    Returns:
        bool: True if M is a point of continuity, False otherwise.
    """
    print("--- Checking Matrix ---")
    print(M)
    
    try:
        if M.shape[0] != M.shape[1]:
            print("Matrix must be square.")
            return False
        n = M.shape[0]
        if n == 0:
            return True # Or handle as an error
            
        # Get eigenvalues. Due to floating point arithmetic, eigenvalues that are
        # theoretically identical might be slightly different. We need to group them.
        eigenvalues = np.linalg.eigvals(M)
        print(f"\nEigenvalues: {eigenvalues}")

        # Group close eigenvalues to find numerically distinct ones
        tol = 1e-8
        distinct_eigenvalues = []
        sorted_evals = np.sort(eigenvalues)
        i = 0
        while i < n:
            lam = sorted_evals[i]
            distinct_eigenvalues.append(lam)
            j = i + 1
            while j < n and np.abs(sorted_evals[j] - lam) < tol:
                j += 1
            i = j
        
        print(f"Distinct eigenvalues: {distinct_eigenvalues}\n")

        # Check the rank condition for each distinct eigenvalue
        identity = np.eye(n)
        is_non_derogatory = True
        for lam in distinct_eigenvalues:
            A = M - lam * identity
            # The rank is the number of non-zero singular values.
            # matrix_rank is robust for floating point matrices.
            r = np.linalg.matrix_rank(A, tol=tol)
            
            print(f"For eigenvalue lambda = {lam:.4f}:")
            print(f"Rank(M - lambda*I) = {r}")
            print(f"Required rank for continuity: n-1 = {n-1}")

            if r != n - 1:
                is_non_derogatory = False
                print("Condition failed: rank != n-1. Matrix is derogatory.\n")
                break
            else:
                print("Condition passed: rank == n-1.\n")
                
        print(f"Conclusion: The matrix is a point of continuity? {is_non_derogatory}")
        return is_non_derogatory

    except Exception as e:
        print(f"An error occurred: {e}")
        return False

if __name__ == '__main__':
    # Example 1: A derogatory matrix (point of discontinuity)
    # The 2x2 identity matrix. It has one eigenvalue 1 with geometric multiplicity 2.
    # rank(I - 1*I) = rank(0) = 0, which is not n-1=1.
    M1 = np.array([[1.0, 0.0], [0.0, 1.0]])
    is_point_of_continuity(M1)

    print("\n" + "="*40 + "\n")
    
    # Example 2: A non-derogatory matrix (point of continuity)
    # A 2x2 Jordan block. It has one eigenvalue 1 with geometric multiplicity 1.
    # rank(M - 1*I) = rank([[0,1],[0,0]]) = 1, which is n-1=1.
    M2 = np.array([[1.0, 1.0], [0.0, 1.0]])
    is_point_of_continuity(M2)

    print("\n" + "="*40 + "\n")
    
    # Example 3: A non-derogatory matrix with distinct eigenvalues (point of continuity)
    # Any matrix with n distinct eigenvalues is non-derogatory.
    M3 = np.array([[1.0, 2.0], [3.0, 4.0]])
    is_point_of_continuity(M3)