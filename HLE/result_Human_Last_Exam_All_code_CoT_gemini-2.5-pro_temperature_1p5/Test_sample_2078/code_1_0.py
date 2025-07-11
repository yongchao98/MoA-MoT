import collections

def solve_cas_probability():
    """
    Solves the probability problem for the described Complex Adaptive System.
    """

    # 1. Define the system's rules and probabilities

    # Subsystem S1: Markov process transition probabilities
    # P_S1[from_state][to_state] = probability
    # We use a defaultdict to return 0 for non-existent transitions.
    P_S1 = {
        'A': collections.defaultdict(float, {'B': 0.3, 'C': 0.5, 'D': 0.2}),
        'B': collections.defaultdict(float, {'A': 0.4, 'C': 0.6}),
        'C': collections.defaultdict(float, {'A': 0.1, 'B': 0.2, 'D': 0.7}),
        'D': collections.defaultdict(float, {'A': 0.8, 'B': 0.2})
    }

    # Coupling S1 -> S2: Deterministic mapping from S1(t) to S2(t+1)
    S1_to_S2 = {'A': 'X', 'B': 'Y', 'C': 'Z', 'D': 'X'}

    # Coupling S2 -> S3: Probabilistic mapping from S2(t) to S3(t+1)
    # P_S3[from_S2_state][to_S3_state] = probability
    P_S3_given_S2 = {
        'X': {'M': 0.7, 'N': 0.3},
        'Y': {'M': 0.4, 'N': 0.6},
        'Z': {'M': 0.2, 'N': 0.8}
    }

    # 2. Define initial and target states
    s1_initial = 'A'
    s1_target = 'D'
    s2_target = 'Z'
    s3_target = 'N'

    # 3. Use logic to find the unique path for the subsystems
    
    # The target state S2(3) = 'Z' implies that S1(2) must be 'C',
    # based on the S1_to_S2 mapping.
    s1_intermediate_2 = 'C'

    # The S1 path over 3 transitions is S1(0) -> S1(1) -> S1(2) -> S1(3).
    # We know S1(0)=A, S1(2)=C, S1(3)=D. We must find S1(1).
    
    total_prob = 0.0
    
    possible_s1_states = ['A', 'B', 'C', 'D']
    
    # Sum probabilities over all possible intermediate paths.
    # In this problem, constraints will prune this to a single path.
    for s1_intermediate_1 in possible_s1_states:
        # Probability of the first transition for S1: S1(0) -> S1(1)
        p_trans1 = P_S1[s1_initial][s1_intermediate_1]

        # Probability of the second transition for S1: S1(1) -> S1(2)
        p_trans2 = P_S1[s1_intermediate_1][s1_intermediate_2]
        
        # Probability of the third transition for S1: S1(2) -> S1(3)
        p_trans3 = P_S1[s1_intermediate_2][s1_target]

        # If any part of the S1 path is impossible, its probability is 0.
        if p_trans1 == 0 or p_trans2 == 0 or p_trans3 == 0:
            continue

        # Determine the state of S2(2), which influences S3(3)
        # S2(2) is determined by S1(1)
        s2_intermediate_2 = S1_to_S2[s1_intermediate_1]

        # Probability of the S3 outcome S3(3) = 'N', given S2(2)
        p_s3_outcome = P_S3_given_S2[s2_intermediate_2][s3_target]
        
        # The probability of this specific path is the product of all independent steps
        path_prob = p_trans1 * p_trans2 * p_trans3 * p_s3_outcome
        
        # Print the breakdown of the calculation for the valid path
        print("Found a valid path for S1: {} -> {} -> {} -> {}".format(s1_initial, s1_intermediate_1, s1_intermediate_2, s1_target))
        print("The probability calculation is based on the following components:")
        print("1. P(S1: {} -> {}): {}".format(s1_initial, s1_intermediate_1, p_trans1))
        print("2. P(S1: {} -> {}): {}".format(s1_intermediate_1, s1_intermediate_2, p_trans2))
        print("3. P(S1: {} -> {}): {}".format(s1_intermediate_2, s1_target, p_trans3))
        print("4. P(S3={} | S2={}): {}".format(s3_target, s2_intermediate_2, p_s3_outcome))
        
        print("\nThe final equation is:")
        print(f"{p_trans1} * {p_trans2} * {p_trans3} * {p_s3_outcome} = {path_prob}")

        total_prob += path_prob

    print("\n---")
    print(f"The total probability of the system being in state ({s1_target}, {s2_target}, {s3_target}) after 3 transitions is: {total_prob}")


solve_cas_probability()
<<<0.0756>>>