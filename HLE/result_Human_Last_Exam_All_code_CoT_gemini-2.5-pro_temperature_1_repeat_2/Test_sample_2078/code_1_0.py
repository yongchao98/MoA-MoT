def solve_cas_probability():
    """
    Calculates the probability of a complex adaptive system reaching a specific state.
    """
    # S1 transition probabilities
    s1_p = {
        'A': {'B': 0.3, 'C': 0.5, 'D': 0.2},
        'B': {'A': 0.4, 'C': 0.6},
        'C': {'A': 0.1, 'B': 0.2, 'D': 0.7},
        'D': {'A': 0.8, 'B': 0.2}
    }

    # S1 -> S2 deterministic coupling map: s2_{t+1} = f(s1_{t})
    s1_s2_map = {'A': 'X', 'B': 'Y', 'C': 'Z', 'D': 'X'}

    # S2 -> S3 probabilistic coupling: P(s3_{t+1} | s2_{t})
    s2_s3_p = {
        'X': {'M': 0.7, 'N': 0.3},
        'Y': {'M': 0.4, 'N': 0.6},
        'Z': {'M': 0.2, 'N': 0.8}
    }

    initial_state = ('A', 'X', 'M')
    target_state = ('D', 'Z', 'N')
    n_steps = 3

    s1_start = initial_state[0]
    s1_target_t3 = target_state[0]
    s2_target_t3 = target_state[1]
    s3_target_t3 = target_state[2]

    # From target state constraints, determine intermediate states
    # To have s2_t3 = 'Z', s1_t2 must be 'C' because s1_s2_map['C'] == 'Z'
    s1_target_t2 = 'C'

    # The only valid S1 path is A -> B -> C -> D
    # Let's confirm this programmatically
    s1_path = []
    total_prob = 0.0
    
    # Step 1: A -> s1_t1
    for s1_t1, p1 in s1_p[s1_start].items():
        # Step 2: s1_t1 -> s1_t2
        if s1_target_t2 in s1_p.get(s1_t1, {}):
            p2 = s1_p[s1_t1][s1_target_t2]
            # Step 3: s1_t2 -> s1_t3
            if s1_target_t3 in s1_p.get(s1_target_t2, {}):
                p3 = s1_p[s1_target_t2][s1_target_t3]
                
                # This is a valid S1 path.
                path = [s1_start, s1_t1, s1_target_t2, s1_target_t3]
                s1_path_prob = p1 * p2 * p3
                
                # Determine s2_t2, which influences the final s3 transition
                s2_t2 = s1_s2_map[s1_t1]
                
                # Get the probability of the final S3 transition
                s3_final_prob = s2_s3_p[s2_t2][s3_target_t3]

                # The total probability for this path
                path_total_prob = s1_path_prob * s3_final_prob
                total_prob += path_total_prob
                
                # Printing the equation for this valid path
                print("Found a valid path of states for S1: {}".format(" -> ".join(path)))
                print("The overall probability is the product of the S1 path probability and the final S3 transition probability.")
                print("P(Total) = P(A->B) * P(B->C) * P(C->D) * P(S3->N | S2 at t=2 was Y)")
                print(f"P(Total) = {p1} * {p2} * {p3} * {s3_final_prob}")
                print(f"P(Total) = {s1_path_prob:.3f} * {s3_final_prob}")
                print(f"P(Total) = {path_total_prob}")


    if total_prob == 0.0:
        print("No valid path found. The probability is 0.")
    else:
        print("\nFinal calculated probability:")
        print(total_prob)

solve_cas_probability()
<<<0.0756>>>