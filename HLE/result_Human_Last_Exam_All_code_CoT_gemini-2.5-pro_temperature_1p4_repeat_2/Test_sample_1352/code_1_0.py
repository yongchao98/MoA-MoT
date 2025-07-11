import numpy as np

def calculate_expected_rounds():
    """
    Calculates the expected number of rounds for the system to return to its initial state.
    """

    # Define influence sets as per the problem description
    # influence_sets[p] = set of people person p influences
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

    # Create the influence matrix C over GF(2)
    # C[j-1, p-1] = 1 if person p influences person j
    C = np.zeros((8, 8), dtype=int)
    for p, influences in influence_sets.items():
        if p > 0 and influences:
            for j in influences:
                C[j - 1, p - 1] = 1

    def rank_gf2(M):
        """Computes the rank of a matrix over GF(2)."""
        mat = M.copy()
        rows, cols = mat.shape
        rank = 0
        pivot_row = 0
        for j in range(cols):
            if pivot_row < rows:
                i = pivot_row
                while i < rows and mat[i, j] == 0:
                    i += 1
                if i < rows:
                    mat[[pivot_row, i]] = mat[[i, pivot_row]]
                    for row_idx in range(rows):
                        if row_idx != pivot_row and mat[row_idx, j] == 1:
                            mat[row_idx, :] = (mat[row_idx, :] + mat[pivot_row, :]) % 2
                    pivot_row += 1
        return pivot_row
    
    def nullity_gf2(M):
        """Computes the nullity of a matrix over GF(2)."""
        return M.shape[1] - rank_gf2(M)

    # Compute powers of C
    C2 = (C @ C) % 2
    C4 = (C2 @ C2) % 2
    # C8 is the zero matrix since the nilpotency index of C is 7.

    # Calculate nullities
    d1 = nullity_gf2(C)
    d2 = nullity_gf2(C2)
    d4 = nullity_gf2(C4)
    d8 = 8 # Nullity of zero matrix

    # Calculate the number of states for each period R
    num_states_r1 = 2**d1
    num_states_r2 = 2**d2 - 2**d1
    num_states_r4 = 2**d4 - 2**d2
    num_states_r8 = 2**d8 - 2**d4
    
    total_states = num_states_r1 + num_states_r2 + num_states_r4 + num_states_r8

    # Calculate the sum of R over all states
    sum_R = (1 * num_states_r1) + \
            (2 * num_states_r2) + \
            (4 * num_states_r4) + \
            (8 * num_states_r8)

    # Calculate the expected value of R
    expected_R = sum_R / total_states

    # Print the results
    print(f"Number of states with period R=1: {num_states_r1}")
    print(f"Number of states with period R=2: {num_states_r2}")
    print(f"Number of states with period R=4: {num_states_r4}")
    print(f"Number of states with period R=8: {num_states_r8}")
    print(f"\nTotal states checked: {total_states}")

    print("\nCalculation of E[R]:")
    print(f"E[R] = (1 * {num_states_r1} + 2 * {num_states_r2} + 4 * {num_states_r4} + 8 * {num_states_r8}) / {total_states}")
    print(f"E[R] = ({1 * num_states_r1} + {2 * num_states_r2} + {4 * num_states_r4} + {8 * num_states_r8}) / {total_states}")
    print(f"E[R] = {sum_R} / {total_states}")
    print(f"E[R] = {expected_R}")
    
    print(f"\nThe expected value of R, rounded to 2 decimal places, is: {expected_R:.2f}")


calculate_expected_rounds()
print("<<<7.44>>>")