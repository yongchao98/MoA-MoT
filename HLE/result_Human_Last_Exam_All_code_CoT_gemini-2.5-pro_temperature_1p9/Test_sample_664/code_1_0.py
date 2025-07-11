def solve_checkerboard_symmetry():
    """
    This function calculates the number of ways to place 8 chips on an 8x8 board
    with one chip per row and column, such that the arrangement is symmetric
    about at least one of the two main diagonals.
    """

    # Let M be the set of configurations symmetric about the main diagonal.
    # |M| is the number of involutions on 8 elements, I_8.
    # The recurrence is I_n = I_{n-1} + (n-1) * I_{n-2}.
    I = {0: 1, 1: 1}
    for n in range(2, 9):
        I[n] = I[n - 1] + (n - 1) * I[n - 2]
    num_main_symmetry = I[8]

    # Let A be the set of configurations symmetric about the anti-diagonal.
    # The recurrence for |A| on an n x n board (n must be even) is
    # A_n = 2*A_{n-2} + (n-2)*A_{n-4}.
    A = {0: 1, 2: 2}
    for n in range(4, 9, 2):
        A[n] = 2 * A[n - 2] + (n - 2) * A[n - 4]
    num_anti_symmetry = A[8]

    # Let |M_intersect_A| be the number of configurations symmetric about both diagonals.
    # As reasoned in the plan, A is a subset of M, so |M_intersect_A| = |A|.
    # We can also compute this directly. A permutation symmetric about both diagonals
    # must be an involution that also commutes with the anti-diagonal reflection.
    # The recurrence for n pairs of elements is U_n = 2*U_{n-1} + 2*(n-1)*U_{n-2}.
    # For an 8x8 board, there are 4 pairs of elements.
    U = {0: 1, 1: 2}
    for n in range(2, 5): # n is the number of pairs
        U[n] = 2 * U[n - 1] + 2 * (n - 1) * U[n - 2]
    num_both_symmetry = U[4]
    
    # The total number of configurations is given by the Principle of Inclusion-Exclusion:
    # Total = |M| + |A| - |M intersect A|
    # Since A is a subset of M, |M intersect A| = |A|. So Total = |M|.
    total_configurations = num_main_symmetry + num_anti_symmetry - num_both_symmetry

    print(f"Number of configurations with main diagonal symmetry: {num_main_symmetry}")
    print(f"Number of configurations with anti-diagonal symmetry: {num_anti_symmetry}")
    print(f"Number of configurations with both symmetries: {num_both_symmetry}")
    print("The total number of possible configurations is calculated by the Principle of Inclusion-Exclusion:")
    print(f"{num_main_symmetry} + {num_anti_symmetry} - {num_both_symmetry} = {total_configurations}")

solve_checkerboard_symmetry()