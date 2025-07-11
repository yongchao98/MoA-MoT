import numpy as np

def check_matrices():
    """
    Checks which of the given matrices are in the set P.

    The set P is the convex hull of matrices v * v.T where v is a non-zero 2D integer vector.
    A matrix M is in P if it satisfies these necessary conditions:
    1. M is symmetric.
    2. M is positive semidefinite.
    3. The trace of M is >= 1.
    And if it can be expressed as a convex combination of the generating matrices.
    """
    # Define the matrices
    A = np.array([[0, 0], [0, 0]])
    B = np.array([[6, 4], [3, 7]])
    C = np.array([[1, -0.5], [-0.5, 1]])
    D = np.array([[np.pi, 1], [1, np.pi**2]])
    E = np.array([[1, np.pi], [np.pi, 1]])
    F = np.array([[42, 0], [0, 0]])
    
    matrices = {'A': A, 'B': B, 'C': C, 'D': D, 'E': E, 'F': F}
    in_P = []

    # Analysis for each matrix
    
    # Matrix A
    # Reason: Fails trace condition (Tr(A) = 0 < 1)
    
    # Matrix B
    # Reason: Not symmetric
    
    # Matrix C
    # It is symmetric, PSD, Tr(C) = 2 >= 1.
    # It can be written as a convex combination:
    # C = (3/4) * v1*v1.T + (1/4) * v2*v2.T
    # where v1 = (1, -1) and v2 = (1, 1).
    v1 = np.array([[1], [-1]])
    v2 = np.array([[1], [1]])
    v1_v1t = v1 @ v1.T
    v2_v2t = v2 @ v2.T
    c_reconstructed = 0.75 * v1_v1t + 0.25 * v2_v2t
    if np.allclose(C, c_reconstructed):
        in_P.append('C')

    # Matrix D
    # Reason: Contains transcendental numbers and cannot be formed from a convex
    # combination of integer matrices. This relies on the fact that pi is transcendental.
    
    # Matrix E
    # Reason: Not positive semidefinite (determinant is 1 - pi^2 < 0)

    # Matrix F
    # It is symmetric, PSD, Tr(F) = 42 >= 1.
    # It can be written as a convex combination:
    # F = (7/13) * v1*v1.T + (6/13) * v2*v2.T
    # where v1 = (6, 0) and v2 = (7, 0).
    v1 = np.array([[6], [0]])
    v2 = np.array([[7], [0]])
    v1_v1t = v1 @ v1.T
    v2_v2t = v2 @ v2.T
    f_reconstructed = (7/13) * v1_v1t + (6/13) * v2_v2t
    if np.allclose(F, f_reconstructed):
        in_P.append('F')
        
    print(in_P)

check_matrices()