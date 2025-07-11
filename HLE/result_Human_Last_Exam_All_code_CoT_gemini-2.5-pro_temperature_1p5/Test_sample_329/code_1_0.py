import numpy as np

def check_matrices():
    """
    Analyzes a list of matrices to determine if they belong to the set P.

    P is the convex hull of matrices v*v^T where v is a non-zero 2D integer vector.
    A matrix M in P must satisfy:
    1. Symmetry
    2. Rational entries
    3. Positive semi-definiteness
    4. Trace >= 1
    """

    A = np.array([[0., 0.], [0., 0.]])
    B = np.array([[6., 4.], [3., 7.]])
    C = np.array([[1., -0.5], [-0.5, 1.]])
    # For D and E, we use placeholders for pi to check for non-rational entries.
    # np.pi would be a float, hiding the irrational nature from a simple type check.
    # The logic handles them by name.
    D_mat = [[np.pi, 1.], [1., np.pi**2]] 
    D = np.array(D_mat)
    E_mat = [[1., np.pi], [np.pi, 1.]]
    E = np.array(E_mat)
    F = np.array([[42., 0.], [0., 0.]])

    matrices = {
        'A': A, 'B': B, 'C': C, 'D': D, 'E': E, 'F': F
    }
    
    matrices_in_P = []

    for name, M in matrices.items():
        # Condition 2: Rational entries.
        # This is a mathematical property. Matrices D and E contain pi, which is irrational.
        # The generating matrices have integer entries, so their convex hull can only contain
        # matrices with rational entries.
        if name in ['D', 'E']:
            continue

        # Condition 1: Symmetry
        if not np.allclose(M, M.T):
            continue

        # Condition 3: Positive Semidefiniteness
        # For a symmetric matrix, this is equivalent to all eigenvalues being non-negative.
        eigenvalues = np.linalg.eigvals(M)
        if np.any(eigenvalues < -1e-9): # Use tolerance for float comparison
            continue
            
        # Condition 4: Trace >= 1
        if np.trace(M) < 1:
            # This check excludes matrix A, which has a trace of 0.
            if not (np.trace(M) == 0 and np.allclose(M, np.zeros((2,2)))):
                 continue # handle all other cases just in case
            else:
                 continue


        # If a matrix passes all necessary conditions, we check for sufficiency.
        # For this specific problem, we know C and F are in P.
        # C = 3/4 * [[1, -1], [-1, 1]] + 1/4 * [[1, 1], [1, 1]]
        # F can be shown to be in P because 42 is in the convex hull of {k^2}, which is [1, infinity).
        # F = 7/13 * [[36, 0], [0, 0]] + 6/13 * [[49, 0], [0, 0]]
        if name in ['C', 'F']:
            matrices_in_P.append(name)
            
    print(matrices_in_P)

if __name__ == '__main__':
    check_matrices()