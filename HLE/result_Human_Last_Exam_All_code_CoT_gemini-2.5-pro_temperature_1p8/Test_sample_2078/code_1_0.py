import collections

def solve():
    """
    Calculates the probability of the system being in state (D, Z, N) after 3 transitions.
    """
    # S1 transition probabilities
    P_S1 = {
        'A': {'B': 0.3, 'C': 0.5, 'D': 0.2},
        'B': {'A': 0.4, 'C': 0.6},
        'C': {'A': 0.1, 'B': 0.2, 'D': 0.7},
        'D': {'A': 0.8, 'B': 0.2}
    }

    # S1 to S2 deterministic mapping
    S1_to_S2 = {
        'A': 'X',
        'B': 'Y',
        'C': 'Z',
        'D': 'X'
    }

    # S3 transition probabilities given S2 state
    P_S3_given_S2 = {
        'X': {'M': 0.7, 'N': 0.3},
        'Y': {'M': 0.4, 'N': 0.6},
        'Z': {'M': 0.2, 'N': 0.8}
    }

    # Initial and target states
    initial_state = ('A', 'X', 'M')
    target_state = ('D', 'Z', 'N')
    
    # Unpack states for clarity
    s1_t0, s2_t0, s3_t0 = initial_state
    s1_t3, s2_t3, s3_t3 = target_state

    # Step 1: Determine the S1 state at t=2.
    # To have s2(t=3) = 'Z', s1(t=2) must be 'C'.
    s1_t2 = 'C'

    # Step 2: Determine the S1 state at t=1.
    # The S1 path must be s1(t=0) -> s1(t=1) -> s1(t=2), which is 'A' -> s1(t=1) -> 'C'.
    # Based on the transition matrix P_S1, the only valid intermediate state is 'B'.
    s1_t1 = 'B'

    # Probabilities for the mandatory S1 path: A -> B -> C -> D
    p_a_b = P_S1[s1_t0][s1_t1]
    p_b_c = P_S1[s1_t1][s1_t2]
    p_c_d = P_S1[s1_t2][s1_t3]

    # Step 3: Determine the S2 state at t=2.
    # This state determines the final S3 transition. It's dictated by s1(t=1).
    s2_t2 = S1_to_S2[s1_t1]

    # Step 4: Determine the probability of the final S3 transition.
    # This is P(S3 becomes N | S2 was Y).
    prob_s3_final_step = P_S3_given_S2[s2_t2][s3_t3]

    # Step 5: Calculate the final probability.
    final_prob = p_a_b * p_b_c * p_c_d * prob_s3_final_step
    
    # Output the explanation and the final equation
    print("The only path for S1 that results in S1(t=3)=D and S2(t=3)=Z is A -> B -> C -> D.")
    print("The S2 state at t=2 is determined by S1 at t=1 ('B'), so S2(t=2)='Y'.")
    print("This determines the probability of the final S3 transition to 'N'.")
    print("\nThe total probability is the product of the probabilities of each required transition:")
    print(f"P = P(A->B) * P(B->C) * P(C->D) * P(S3=N|S2=Y)")
    print(f"P = {p_a_b} * {p_b_c} * {p_c_d} * {prob_s3_final_step} = {final_prob:.4f}")
    
    return final_prob

# Run the solver
solve()