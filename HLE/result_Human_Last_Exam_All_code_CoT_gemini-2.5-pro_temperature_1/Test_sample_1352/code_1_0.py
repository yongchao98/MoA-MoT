import numpy as np

def solve():
    """
    Calculates the expected number of rounds for the system to return to its initial state.
    """
    # Step 1: Define influence sets and create the influence matrix M.
    # M[i, j] = 1 if person j+1 influences person i+1.
    influence = {
        1: {2, 4, 6, 7}, 2: {3, 5, 6, 8}, 3: {4, 6},   4: {5},
        5: {6, 8},      6: {7},           7: {8},      8: {}
    }
    M = np.zeros((8, 8), dtype=int)
    for influencer_idx in range(1, 9):
        for influenced_idx in influence[influencer_idx]:
            # Matrix indices are 0-based
            M[influenced_idx - 1, influencer_idx - 1] = 1

    # Helper function to compute matrix rank over F_2 (mod 2).
    def rank_mod2(A):
        m, n = A.shape
        mat = A.copy().astype(int)
        rank = 0
        pivot_row = 0
        for j in range(n):  # For each column
            if pivot_row < m:
                pivot = pivot_row
                while pivot < m and mat[pivot, j] == 0:
                    pivot += 1
                if pivot < m:
                    mat[[pivot_row, pivot]] = mat[[pivot, pivot_row]]
                    for i in range(m):
                        if i != pivot_row and mat[i, j] == 1:
                            mat[i, :] = (mat[i, :] + mat[pivot_row, :]) % 2
                    pivot_row += 1
        return pivot_row

    # Step 2: Compute the required powers of M. The order of T is 8,
    # so we need to check orders 1, 2, 4, 8. This requires M, M^2, M^4, M^8.
    M2 = (M @ M) % 2
    M4 = (M2 @ M2) % 2
    M8 = (M4 @ M4) % 2

    # Step 3: Compute the ranks of these matrices.
    rank_M1 = rank_mod2(M)
    rank_M2 = rank_mod2(M2)
    rank_M4 = rank_mod2(M4)
    rank_M8 = rank_mod2(M8) # This should be 0, as M is nilpotent of index 6.

    # Step 4: Calculate the size of the kernel for each matrix.
    # |ker(A)| = 2^(n - rank(A)) where n=8.
    ker_M1_size = 2**(8 - rank_M1)
    ker_M2_size = 2**(8 - rank_M2)
    ker_M4_size = 2**(8 - rank_M4)
    ker_M8_size = 2**(8 - rank_M8)

    # Step 5: Count the number of states for each possible order R.
    # C_d is the number of states with order exactly d.
    C1 = ker_M1_size
    C2 = ker_M2_size - ker_M1_size
    C4 = ker_M4_size - ker_M2_size
    C8 = ker_M8_size - ker_M4_size

    # Step 6: Calculate the sum of all R values over the 256 states.
    total_R = (1 * C1) + (2 * C2) + (4 * C4) + (8 * C8)
    
    # Step 7: Calculate the expected value E[R].
    expected_R = total_R / 256.0

    # Step 8: Print the breakdown of the calculation and the final answer.
    print("The number of initial states returning after exactly d rounds:")
    print(f"d=1: {C1} states")
    print(f"d=2: {C2} states")
    print(f"d=4: {C4} states")
    print(f"d=8: {C8} states")
    print("\nThe expected value E[R] is the average of R over all 256 states:")
    print(f"E[R] = (1 * {C1} + 2 * {C2} + 4 * {C4} + 8 * {C8}) / 256")
    print(f"E[R] = ({1 * C1} + {2 * C2} + {4 * C4} + {8 * C8}) / 256")
    print(f"E[R] = {total_R} / 256")
    print(f"E[R] = {expected_R:.2f}")

solve()
<<<6.92>>>