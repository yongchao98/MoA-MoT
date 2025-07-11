import collections

def solve_system_probability():
    """
    Calculates the probability of the system reaching state (D, Z, N) from (A, X, M) in 3 transitions.
    """

    # 1. Define the system's rules and probabilities

    # S1: Markov process transition probabilities
    s1_probs = {
        'A': {'B': 0.3, 'C': 0.5, 'D': 0.2},
        'B': {'A': 0.4, 'C': 0.6},
        'C': {'A': 0.1, 'B': 0.2, 'D': 0.7},
        'D': {'A': 0.8, 'B': 0.2},
    }

    # S1-S2 Coupling: s2(t+1) is determined by s1(t)
    s1_s2_map = {'A': 'X', 'B': 'Y', 'C': 'Z', 'D': 'X'}

    # S2-S3 Coupling: s3(t+1) probability is determined by s2(t)
    s2_s3_probs = {
        'X': {'M': 0.7, 'N': 0.3},
        'Y': {'M': 0.4, 'N': 0.6},
        'Z': {'M': 0.2, 'N': 0.8},
    }
    
    # 2. Identify the unique path for S1 and S2
    # Initial state at t=0: S1=A, S2=X, S3=M
    # Target state at t=3: S1=D, S2=Z, S3=N
    
    # Constraint: For s2_3 to be 'Z', s1_2 must be 'C' (since s1_s2_map['C'] == 'Z')
    # So the S1 path must be: s1_0 -> s1_1 -> 'C' -> 'D'
    # From s1_0 = 'A', the only state that can transition to 'C' is 'B'
    # (Checking: D cannot transition to C, C->C prob is 0).
    # Therefore, the unique S1 path is A -> B -> C -> D
    s1_path = ['A', 'B', 'C', 'D']
    
    # This determines the S2 path based on s2(t+1) = map(s1(t))
    # s2_0 = 'X' (initial)
    # s2_1 = map(s1_0='A') = 'X'
    # s2_2 = map(s1_1='B') = 'Y'
    # s2_3 = map(s1_2='C') = 'Z'
    s2_path = ['X', 'X', 'Y', 'Z']
    
    # 3. Calculate the total probability

    # Probability of the S1 path A -> B -> C -> D
    p_b_given_a = s1_probs['A']['B']
    p_c_given_b = s1_probs['B']['C']
    p_d_given_c = s1_probs['C']['D']
    s1_path_prob = p_b_given_a * p_c_given_b * p_d_given_c
    
    # Probability of S3 ending in 'N'.
    # This is determined by the state of S2 at t=2, which is s2_2='Y'.
    # The intermediate states of S3 do not affect this final transition's probability
    # because the probabilities of reaching those intermediate states sum to 1.
    s2_driver_state_for_s3_final = s2_path[2] # state at t=2
    s3_final_prob = s2_s3_probs[s2_driver_state_for_s3_final]['N']
    
    # The total probability is the product of the S1 path probability and the
    # conditional probability of S3 ending in N.
    total_prob = s1_path_prob * s3_final_prob
    
    # Output the explanation and final equation
    print("Step 1: Determine the unique path for Subsystem S1.")
    print(f"The final state requires S1(3)=D and S2(3)=Z.")
    print(f"The S1-S2 coupling (s2(t+1) = map(s1(t))) means S1(2) must be 'C' for S2(3) to be 'Z'.")
    print(f"The only 3-step path for S1 from A to D with the second state being C is {s1_path[0]} -> {s1_path[1]} -> {s1_path[2]} -> {s1_path[3]}.\n")
    
    print("Step 2: Calculate the probability of this S1 path.")
    print(f"P(S1 Path) = P({s1_path[1]}|{s1_path[0]}) * P({s1_path[2]}|{s1_path[1]}) * P({s1_path[3]}|{s1_path[2]})")
    print(f"P(S1 Path) = {p_b_given_a} * {p_c_given_b} * {p_d_given_c} = {s1_path_prob:.3f}\n")
    
    print("Step 3: Determine the required S3 transition probability.")
    print(f"The S1 path determines the S2 path: {s2_path[0]} -> {s2_path[1]} -> {s2_path[2]} -> {s2_path[3]}.")
    print(f"The final transition of S3 (from t=2 to t=3) is driven by the state of S2 at t=2, which is '{s2_driver_state_for_s3_final}'.")
    print(f"The probability of S3 transitioning to 'N' given S2 is in state '{s2_driver_state_for_s3_final}' is {s3_final_prob}.\n")
    
    print("Step 4: Calculate the final probability.")
    print("Total Probability = P(S1 Path) * P(S3(3)='N' | S2(2)='Y')")
    print(f"Total Probability = {s1_path_prob:.3f} * {s3_final_prob}")
    print(f"Total Probability = {total_prob}")
    print("\nThe final equation is:")
    print(f"({p_b_given_a} * {p_c_given_b} * {p_d_given_c}) * {s3_final_prob} = {total_prob}")


solve_system_probability()