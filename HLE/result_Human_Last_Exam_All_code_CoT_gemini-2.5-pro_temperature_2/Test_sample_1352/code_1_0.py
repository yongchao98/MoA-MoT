import numpy as np

def solve():
    """
    Calculates the expected number of rounds for the switches to return to their initial state.
    """
    # Step 1: Define the influence matrix M.
    # M[i, j] = 1 means person j+1 influences person i+1.
    influence_sets = [
        {2, 4, 6, 7},  # Person 1 influences
        {3, 5, 6, 8},  # Person 2 influences
        {4, 6},        # Person 3 influences
        {5},           # Person 4 influences
        {6, 8},        # Person 5 influences
        {7},           # Person 6 influences
        {8},           # Person 7 influences
        {}             # Person 8 influences
    ]
    M = np.zeros((8, 8), dtype=int)
    for influencer_idx, influenced_set in enumerate(influence_sets):
        for influenced_person_num in influenced_set:
            M[influenced_person_num - 1, influencer_idx] = 1

    # The transformation is s_new = (I + M) @ s (mod 2)
    # The cycle length r for a state s satisfies (A^r - I)s = 0, where A = I + M.
    # For r=2^k, this simplifies to M^r @ s = 0.
    # The possible cycle lengths are divisors of the order of A, which is 8.
    
    # Step 2: Compute powers of M and their ranks to find dimensions of kernels.
    # dim(ker(M^d)) = 8 - rank(M^d)
    M1 = M
    M2 = (M @ M) % 2
    M4 = (M2 @ M2) % 2
    # M8 is the zero matrix since M is nilpotent of index 8.

    rank_M1 = np.linalg.matrix_rank(M1)
    rank_M2 = np.linalg.matrix_rank(M2)
    rank_M4 = np.linalg.matrix_rank(M4)
    rank_M8 = 0

    # Step 3: Count states for each cycle length.
    # num_states_div_d is the number of states whose cycle length divides d.
    num_states_div_1 = 2**(8 - rank_M1)
    num_states_div_2 = 2**(8 - rank_M2)
    num_states_div_4 = 2**(8 - rank_M4)
    num_states_div_8 = 2**(8 - rank_M8)

    # N_d is the number of states with cycle length exactly d.
    N_1 = num_states_div_1
    N_2 = num_states_div_2 - num_states_div_1
    N_4 = num_states_div_4 - num_states_div_2
    N_8 = num_states_div_8 - num_states_div_4

    # Step 4: Calculate the expected value E[R].
    total_R_sum = (N_1 * 1) + (N_2 * 2) + (N_4 * 4) + (N_8 * 8)
    expected_R = total_R_sum / 256
    
    print("The expected number of rounds is calculated as follows:")
    print(f"E[R] = (N_1*1 + N_2*2 + N_4*4 + N_8*8) / 256")
    print(f"E[R] = ({N_1} * 1 + {N_2} * 2 + {N_4} * 4 + {N_8} * 8) / 256")
    print(f"E[R] = ({N_1*1} + {N_2*2} + {N_4*4} + {N_8*8}) / 256")
    print(f"E[R] = {total_R_sum} / 256")
    print(f"E[R] = {expected_R:.2f}")

solve()