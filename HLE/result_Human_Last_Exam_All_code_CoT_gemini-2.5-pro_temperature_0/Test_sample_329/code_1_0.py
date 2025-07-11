import numpy as np

def solve_matrices():
    """
    This script analyzes a list of matrices to determine if they belong to the set P,
    which is the convex hull of matrices v*v^T for non-zero 2D integer vectors v.
    It verifies the matrices that are in P by constructing them as a convex combination
    and provides brief reasons for the exclusion of the others.
    """
    result_list = []

    print("--- Analysis of Matrix C ---")
    # C = [[1, -1/2], [-1/2, 1]]
    # C satisfies the necessary properties (Symmetric, PSD, Trace=2).
    # It can be constructed with v1=(1,1) and v2=(1,-1).
    l1_C = 0.25
    S1_C = np.array([[1, 1], [1, 1]])

    l2_C = 0.75
    S2_C = np.array([[1, -1], [-1, 1]])

    C_reconstructed = l1_C * S1_C + l2_C * S2_C
    result_list.append('C')
    print("Matrix C is in P.")
    print("The convex combination is:")
    print(f"{l1_C} * [[{S1_C[0,0]}, {S1_C[0,1]}], [{S1_C[1,0]}, {S1_C[1,1]}]] + {l2_C} * [[{S2_C[0,0]}, {S2_C[0,1]}], [{S2_C[1,0]}, {S2_C[1,1]}]] = [[{C_reconstructed[0,0]}, {C_reconstructed[0,1]}], [{C_reconstructed[1,0]}, {C_reconstructed[1,1]}]]")
    print("-" * 20)


    print("--- Analysis of Matrix F ---")
    # F = [[42, 0], [0, 0]]
    # F satisfies the necessary properties (Symmetric, PSD, Trace=42).
    # The zero entries imply generating vectors must be of the form (x, 0).
    # It can be constructed with v1=(6,0) and v2=(7,0).
    l1_F = 7/13
    S1_F = np.array([[36, 0], [0, 0]])

    l2_F = 6/13
    S2_F = np.array([[49, 0], [0, 0]])

    F_reconstructed = l1_F * S1_F + l2_F * S2_F
    result_list.append('F')
    print("Matrix F is in P.")
    print("The convex combination is:")
    print(f"{l1_F} * [[{S1_F[0,0]}, {S1_F[0,1]}], [{S1_F[1,0]}, {S1_F[1,1]}]] + {l2_F} * [[{S2_F[0,0]}, {S2_F[0,1]}], [{S2_F[1,0]}, {S2_F[1,1]}]] = [[{F_reconstructed[0,0]:.0f}, {F_reconstructed[0,1]:.0f}], [{F_reconstructed[1,0]:.0f}, {F_reconstructed[1,1]:.0f}]]")
    print("-" * 20)

    print("Matrices A, B, D, E are not in P for the following reasons:")
    print("A: Trace is 0, but must be >= 1.")
    print("B: Not symmetric.")
    print("D: Irrational entries (pi) lead to a contradiction with its transcendentality.")
    print("E: Not positive semidefinite (determinant is negative).")
    print("-" * 20)

    print("The final list of matrices in P is:")
    print(result_list)

solve_matrices()