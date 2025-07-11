import numpy as np

def solve_switch_problem():
    """
    Solves the switch problem by modeling it as a linear system over F_2.
    """
    
    # Step 1: Define influence sets and build the action matrix A
    # Person j's influence set I_j is represented by column j-1 of A.
    # A[i, j] = 1 if person j+1 influences person i+1.
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
    
    n = 8
    A = np.zeros((n, n), dtype=int)
    for person_idx, influenced_set in influence_sets.items():
        j = person_idx - 1
        for influenced_person_idx in influenced_set:
            i = influenced_person_idx - 1
            A[i, j] = 1

    # Step 2: Implement a function to calculate matrix rank over F_2
    def rank_F2(M):
        """Calculates the rank of a matrix over the field F_2 using Gaussian elimination."""
        mat = M.copy()
        rows, cols = mat.shape
        rank = 0
        pivot_row = 0
        for j in range(cols):  # Iterate through columns
            if pivot_row < rows:
                i = pivot_row
                while i < rows and mat[i, j] == 0:
                    i += 1
                
                if i < rows:  # Found a pivot at (i, j)
                    # Swap rows to bring pivot to the current pivot_row
                    mat[[pivot_row, i]] = mat[[i, pivot_row]]
                    
                    # Eliminate other 1s in this column
                    for k in range(rows):
                        if k != pivot_row and mat[k, j] == 1:
                            mat[k, :] = (mat[k, :] + mat[pivot_row, :]) % 2
                    
                    pivot_row += 1
        return pivot_row

    # Step 3: Compute powers of A and their ranks to find period distributions
    # The order of matrix M = I+A is 8, so periods can be 1, 2, 4, 8.
    # Number of states with period dividing k=2^p is |Null(A^k)| = 2^(n - rank(A^k))
    A1 = A
    A2 = (A1 @ A) % 2
    A4 = (A2 @ A2) % 2
    A8 = (A4 @ A4) % 2 # This will be the zero matrix

    rank1 = rank_F2(A1)
    rank2 = rank_F2(A2)
    rank4 = rank_F2(A4)
    rank8 = rank_F2(A8) # This will be 0

    num_states_div_1 = 2**(n - rank1)
    num_states_div_2 = 2**(n - rank2)
    num_states_div_4 = 2**(n - rank4)
    num_states_div_8 = 2**(n - rank8)

    # Step 4: Calculate the number of states for each exact period
    N1 = num_states_div_1
    N2 = num_states_div_2 - num_states_div_1
    N4 = num_states_div_4 - num_states_div_2
    N8 = num_states_div_8 - num_states_div_4

    # Step 5: Calculate the expected value E[R] and print the detailed calculation
    total_states = 2**n
    term1 = N1 * 1
    term2 = N2 * 2
    term3 = N4 * 4
    term4 = N8 * 8
    sum_of_periods = term1 + term2 + term3 + term4
    expected_R = sum_of_periods / total_states

    print("The state evolution is a linear transformation. The order of this transformation matrix is 8, so the possible periods for any state are 1, 2, 4, or 8.")
    print("\nCalculating the number of states for each period:")
    print(f"Number of states with period 1 (N_1): {N1}")
    print(f"Number of states with period 2 (N_2): {N2}")
    print(f"Number of states with period 4 (N_4): {N4}")
    print(f"Number of states with period 8 (N_8): {N8}")
    
    print("\nThe expected value E[R] is the average of the periods over all 256 possible initial states.")
    print("E[R] = (1/Total States) * (N_1*1 + N_2*2 + N_4*4 + N_8*8)")
    print(f"E[R] = (1/{total_states}) * ({N1}*1 + {N2}*2 + {N4}*4 + {N8}*8)")
    print(f"E[R] = (1/{total_states}) * ({term1} + {term2} + {term3} + {term4})")
    print(f"E[R] = (1/{total_states}) * ({sum_of_periods})")
    print(f"E[R] = {sum_of_periods / total_states}")
    print(f"\nRounded to 2 decimal places, the final answer is:")
    print(f"{expected_R:.2f}")

solve_switch_problem()
<<<7.46>>>