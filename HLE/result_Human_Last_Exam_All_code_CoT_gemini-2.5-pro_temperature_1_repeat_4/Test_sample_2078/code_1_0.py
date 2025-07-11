import sys

def solve():
    """
    Calculates the probability of a complex system being in a specific state after 3 transitions.
    """
    # Define system parameters based on the problem description
    p_s1 = {
        'A': {'B': 0.3, 'C': 0.5, 'D': 0.2},
        'B': {'A': 0.4, 'C': 0.6},
        'C': {'A': 0.1, 'B': 0.2, 'D': 0.7},
        'D': {'A': 0.8, 'B': 0.2}
    }
    s1_to_s2 = {'A': 'X', 'B': 'Y', 'C': 'Z', 'D': 'X'}
    p_s3_given_s2 = {
        'X': {'M': 0.7, 'N': 0.3},
        'Y': {'M': 0.4, 'N': 0.6},
        'Z': {'M': 0.2, 'N': 0.8}
    }
    
    # Initial and Target states
    s1_initial = 'A'
    s1_target, s2_target, s3_target = 'D', 'Z', 'N'

    # Based on the plan, the required S1 path is A -> B -> C -> D.
    # We extract the probabilities for each step in this path.
    p_A_to_B = p_s1.get(s1_initial, {}).get('B', 0)
    p_B_to_C = p_s1.get('B', {}).get('C', 0)
    p_C_to_D = p_s1.get('C', {}).get(s1_target, 0)
    
    # The probability of this unique S1 path (A -> B -> C -> D)
    prob_s1_path = p_A_to_B * p_B_to_C * p_C_to_D
    
    # The final S3 transition probability depends on S2 at t=2.
    # S2_2 is determined by S1_1. In our path, S1_1 is 'B', so S2_2 must be 'Y'.
    s1_state_at_t1 = 'B'
    s2_state_at_t2 = s1_to_s2[s1_state_at_t1]
    
    # The probability of S3 transitioning to 'N' given S2 was 'Y'
    prob_s3_trans = p_s3_given_s2[s2_state_at_t2][s3_target]

    # The final probability is the product of the S1 path probability and the S3 transition probability.
    total_prob = prob_s1_path * prob_s3_trans
    
    print("The final probability is the product of the probability of the required S1 path and the probability of the required final S3 transition.")
    print("P_total = P(S1 path: A->B->C->D) * P(S3_3=N | S2_2=Y)")
    print(f"P_total = (P(A->B) * P(B->C) * P(C->D)) * P(S3_3=N | S2_2=Y)")
    print(f"P_total = ({p_A_to_B} * {p_B_to_C} * {p_C_to_D}) * {prob_s3_trans}")
    print(f"P_total = {prob_s1_path:.3f} * {prob_s3_trans}")
    print(f"P_total = {total_prob}")
    
    print(f"<<<{total_prob}>>>")

solve()