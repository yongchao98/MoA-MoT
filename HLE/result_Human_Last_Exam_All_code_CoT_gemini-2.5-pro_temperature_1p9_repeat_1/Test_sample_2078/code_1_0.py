import collections

def solve_cas_probability():
    """
    Calculates the probability of a complex adaptive system reaching a specific state
    after three transitions.
    """
    # S1: Transition probabilities
    P_S1 = {
        'A': {'B': 0.3, 'C': 0.5, 'D': 0.2},
        'B': {'A': 0.4, 'C': 0.6},
        'C': {'A': 0.1, 'B': 0.2, 'D': 0.7},
        'D': {'A': 0.8, 'B': 0.2}
    }

    # S1 -> S2 deterministic coupling: s2(t+1) = f12(s1(t))
    f12 = {'A': 'X', 'B': 'Y', 'C': 'Z', 'D': 'X'}

    # S2 -> S3 probabilistic coupling: P(s3(t+1) | s2(t))
    g23 = {
        'X': {'M': 0.7, 'N': 0.3},
        'Y': {'M': 0.4, 'N': 0.6},
        'Z': {'M': 0.2, 'N': 0.8}
    }
    
    # Initial state (t=0)
    s0 = ('A', 'X', 'M')
    # Target state (t=3)
    target_s3 = ('D', 'Z', 'N')

    # From backward reasoning, we know the path for s1 must be A -> B -> C -> D.
    # The path for s2 is then determined: ? -> X -> Y -> Z.
    # This leaves four possible paths depending on the intermediate states of s3.

    # Paths are defined by the state sequence: s0 -> s1 -> s2 -> s3
    # where si = (s1_i, s2_i, s3_i)
    paths_s3_sequences = [('M', 'M'), ('M', 'N'), ('N', 'M'), ('N', 'N')]

    path_probabilities = []
    
    # Calculate transition probabilities for T1, T2, T3
    # P(s(t+1)|s(t)) = P(s1(t)->s1(t+1)) * P(s3(t+1)|s2(t))
    
    # Path states
    s1_path = ['A', 'B', 'C', 'D']
    s2_path = [f12[s1_path[0]], f12[s1_path[1]], f12[s1_path[2]]]
    # Initial state s2(0) is X
    s2_t_minus_1_path = ['X', s2_path[0], s2_path[1]]

    prob_details = []

    for s3_1_state, s3_2_state in paths_s3_sequences:
        s3_t_minus_1_path = [s0[2], s3_1_state, s3_2_state]

        # Transition t=0 to t=1
        # from (A, X, M) to (B, X, s3_1_state)
        p_t1 = P_S1[s1_path[0]][s1_path[1]] * g23[s2_t_minus_1_path[0]][s3_1_state]
        
        # Transition t=1 to t=2
        # from (B, X, s3_1_state) to (C, Y, s3_2_state)
        p_t2 = P_S1[s1_path[1]][s1_path[2]] * g23[s2_t_minus_1_path[1]][s3_2_state]
        
        # Transition t=2 to t=3
        # from (C, Y, s3_2_state) to (D, Z, N)
        p_t3 = P_S1[s1_path[2]][s1_path[3]] * g23[s2_t_minus_1_path[2]][target_s3[2]]

        path_prob = p_t1 * p_t2 * p_t3
        path_probabilities.append(path_prob)
        prob_details.append({
            'path_desc': f"(A,X,M)->(B,X,{s3_1_state})->(C,Y,{s3_2_state})->(D,Z,N)",
            'p_t1': p_t1,
            'p_t2': p_t2,
            'p_t3': p_t3,
            'total': path_prob
        })

    total_prob = sum(path_probabilities)

    print("The final probability is the sum of probabilities of all valid paths from (A,X,M) to (D,Z,N).\n")
    print("There are 4 such paths, differing only in the state of subsystem S3.")
    print("The probability of a path S(0)->S(1)->S(2)->S(3) is P(S(1)|S(0)) * P(S(2)|S(1)) * P(S(3)|S(2)).\n")

    path_eq_parts = []
    for i, detail in enumerate(prob_details):
        print(f"Path {i+1}: {detail['path_desc']}")
        # Using :.4g to format numbers for readability without excessive decimal places.
        print(f"  P(S(1)|S(0)) = P(A->B) * P(S3= {detail['path_desc'][15]} |S2=X) = {P_S1['A']['B']} * {g23['X'][detail['path_desc'][15]]} = {detail['p_t1']:.4g}")
        print(f"  P(S(2)|S(1)) = P(B->C) * P(S3= {detail['path_desc'][22]} |S2=X) = {P_S1['B']['C']} * {g23['X'][detail['path_desc'][22]]} = {detail['p_t2']:.4g}")
        print(f"  P(S(3)|S(2)) = P(C->D) * P(S3=N|S2=Y) = {P_S1['C']['D']} * {g23['Y']['N']} = {detail['p_t3']:.4g}")
        print(f"  Path {i+1} Probability = {detail['p_t1']:.4g} * {detail['p_t2']:.4g} * {detail['p_t3']:.4g} = {detail['total']:.6f}\n")
        path_eq_parts.append(f"{detail['total']:.6f}")

    equation = " + ".join(path_eq_parts)
    print("Final Calculation:")
    print(f"Total Probability = {equation}")
    print(f"Total Probability = {total_prob:.4f}")

solve_cas_probability()
<<<0.0756>>>