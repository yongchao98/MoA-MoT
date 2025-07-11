import collections

def solve_system_probability():
    """
    Calculates the probability of the system being in state (D, Z, N) after 3 transitions.
    """
    # --- Define System Rules ---

    # S1: Markov process transition probabilities
    p_s1 = {
        'A': {'B': 0.3, 'C': 0.5, 'D': 0.2},
        'B': {'A': 0.4, 'C': 0.6},
        'C': {'A': 0.1, 'C': 0.0, 'B': 0.2, 'D': 0.7}, # C->C is not given, assume 0
        'D': {'A': 0.8, 'B': 0.2}
    }

    # S1 -> S2 coupling: S2's state at t+1 is determined by S1's state at t.
    s1_to_s2 = {
        'A': 'X',
        'B': 'Y',
        'C': 'Z',
        'D': 'X'  # Assumed completion of the "and so on" rule
    }

    # S2 -> S3 coupling: S3's state at t+1 is determined by S2's state at t.
    p_s3_given_s2 = {
        'X': {'M': 0.7, 'N': 0.3},
        'Y': {'M': 0.4, 'N': 0.6},
        'Z': {'M': 0.2, 'N': 0.8}
    }

    # --- Initialize State ---
    # probs stores { (s1, s2, s3): probability }
    probs = collections.defaultdict(float)
    initial_state = ('A', 'X', 'M')
    probs[initial_state] = 1.0

    # --- Run Transitions ---
    for t in range(3):
        next_probs = collections.defaultdict(float)
        for state, prob in probs.items():
            if prob == 0:
                continue
            
            s1_t, s2_t, s3_t = state

            # S1 transitions
            for s1_tp1, p_s1_trans in p_s1[s1_t].items():
                # S2's next state is determined by s1_t
                s2_tp1 = s1_to_s2[s1_t]
                
                # S3's next state is determined by s2_t
                for s3_tp1, p_s3_trans in p_s3_given_s2[s2_t].items():
                    # Calculate probability of this specific transition path
                    transition_prob = p_s1_trans * p_s3_trans
                    
                    # Add to the probability of the resulting state
                    new_state = (s1_tp1, s2_tp1, s3_tp1)
                    next_probs[new_state] += prob * transition_prob
        
        probs = next_probs
        # At t=1 (first iteration, t=0), we have probabilities for t=1 states.
        # We need probabilities at t=2 to calculate the final step.
        if t == 1:
            prob_C_Y_M_t2 = probs[('C', 'Y', 'M')]
            prob_C_Y_N_t2 = probs[('C', 'Y', 'N')]

    # --- Calculate and Print Final Answer ---
    target_state = ('D', 'Z', 'N')
    
    # The target state (D, Z, N) at t=3 requires S1 was C at t=2.
    # The only states at t=2 with S1='C' that have non-zero probability are ('C', 'Y', 'M') and ('C', 'Y', 'N').
    
    # Contribution from state ('C', 'Y', 'M') at t=2
    # Transition: C -> D for S1
    # S2 state at t=3 is determined by S1 at t=2 (C), so S2 becomes Z.
    # S3 state at t=3 is determined by S2 at t=2 (Y), so S3 becomes N.
    p_trans_from_Y = p_s1['C']['D'] * p_s3_given_s2['Y']['N']
    contrib1 = prob_C_Y_M_t2 * p_trans_from_Y

    # Contribution from state ('C', 'Y', 'N') at t=2
    contrib2 = prob_C_Y_N_t2 * p_trans_from_Y

    final_prob = probs[target_state]

    print("The final state (D, Z, N) can only be reached from states at t=2 where S1 is 'C'.")
    print("The relevant states at t=2 are ('C', 'Y', 'M') and ('C', 'Y', 'N').")
    print(f"\nProbability of being in state ('C', 'Y', 'M') at t=2 is: {prob_C_Y_M_t2:.4f}")
    print(f"Probability of being in state ('C', 'Y', 'N') at t=2 is: {prob_C_Y_N_t2:.4f}")
    
    p_C_to_D = p_s1['C']['D']
    p_N_given_Y = p_s3_given_s2['Y']['N']
    
    print("\nThe transition probability from these states to the target state involves:")
    print(f"1. S1 transitioning from C to D (Prob: {p_C_to_D})")
    print(f"2. S3 transitioning to N, given S2 was Y (Prob: {p_N_given_Y})")
    
    print("\nThe final probability is the sum of the probabilities of these two paths:")
    print(f"Path 1: P(C,Y,M at t=2) * P(C->D) * P(N|Y) = {prob_C_Y_M_t2:.4f} * {p_C_to_D} * {p_N_given_Y} = {contrib1:.5f}")
    print(f"Path 2: P(C,Y,N at t=2) * P(C->D) * P(N|Y) = {prob_C_Y_N_t2:.4f} * {p_C_to_D} * {p_N_given_Y} = {contrib2:.5f}")
    print(f"\nTotal Probability = {contrib1:.5f} + {contrib2:.5f}")
    print(f"The probability that the system will be in the state (D, Z, N) after 3 transitions is: {final_prob:.4f}")

solve_system_probability()