#
# Name: complex_system_probability.py
#
# Description: Solves the probability of a complex adaptive system reaching a specific state.
#
def solve_system_probability():
    """
    Calculates the probability of a complex adaptive system reaching a target state
    by analyzing the state transition paths based on the problem's rules.
    """

    # Define S1 transition probabilities
    s1_probs = {
        'A': {'B': 0.3, 'C': 0.5, 'D': 0.2},
        'B': {'A': 0.4, 'C': 0.6},
        'C': {'A': 0.1, 'B': 0.2, 'D': 0.7},
        'D': {'A': 0.8, 'B': 0.2}
    }

    # Define S2's state mapping from S1's previous state, s2(t+1) = g(s1(t))
    # Assuming the "and so on" implies a mapping: A->X, B->Y, C->Z, and D wraps back to X.
    s2_mapping = {'A': 'X', 'B': 'Y', 'C': 'Z', 'D': 'X'}

    # Define S3 transition probabilities based on S2's previous state
    # P(s3(t+1) = state | s2(t))
    s3_probs = {
        'X': {'M': 0.7, 'N': 0.3},
        'Y': {'M': 0.4, 'N': 0.6},
        'Z': {'M': 0.2, 'N': 0.8}
    }

    # Initial and Target States
    s1_initial, s2_initial, s3_initial = ('A', 'X', 'M')
    s1_target, s2_target, s3_target = ('D', 'Z', 'N')

    print("This problem asks for the probability of reaching state (D, Z, N) after 3 transitions from (A, X, M).")
    print("The solution involves deducing the single valid path the system must follow to reach the target.")
    
    # --- Step 1: Determine the required path for S1 ---
    print("\nStep 1: Determine the required path for subsystem S1.")
    print("The system evolves from t=0 to t=3.")
    print(f"The target s2(3)='{s2_target}' is determined by s1(2).")
    s1_at_t2 = 'C' # s2_mapping['C'] == 'Z'
    print(f"For s2(3) to be 'Z', s1(2) must be 'C' based on the mapping rule.")

    print(f"\nNow we trace back. For s1(2) to be 'C', s1(1) must be a state that can transition to 'C'.")
    print("The initial state is s1(0)='A'. To get to s1(2)='C', the path s1(0)->s1(1)->s1(2) must be 'A'->s1(1)->'C'.")
    s1_at_t1 = 'B' # P(A->A)=0, so s1(1) can't be A. P(A->B) exists, and P(B->C) exists. This is the only way.
    print("The only valid intermediate state s1(1) is 'B'.")
    
    print(f"\nThe full required path for S1 is: {s1_initial} -> {s1_at_t1} -> {s1_at_t2} -> {s1_target}.")
    p_s1_step1 = s1_probs[s1_initial][s1_at_t1]
    p_s1_step2 = s1_probs[s1_at_t1][s1_at_t2]
    p_s1_step3 = s1_probs[s1_at_t2][s1_target]
    prob_s1_path = p_s1_step1 * p_s1_step2 * p_s1_step3
    print(f"The probability of this S1 path is P({s1_initial}->{s1_at_t1}) * P({s1_at_t1}->{s1_at_t2}) * P({s1_at_t2}->{s1_target}) = {p_s1_step1} * {p_s1_step2} * {p_s1_step3} = {prob_s1_path:.4f}.")

    # --- Step 2: Determine the probability for S3's final transition ---
    print("\nStep 2: Determine the probability of the required S3 transition.")
    print(f"The target s3(3)='{s3_target}' depends on s2(2).")
    s2_at_t2 = s2_mapping[s1_at_t1]
    print(f"The state s2(2) is determined by s1(1). Since the path requires s1(1)='{s1_at_t1}', then s2(2) must be '{s2_at_t2}'.")
    prob_s3_final_step = s3_probs[s2_at_t2][s3_target]
    print(f"The probability of the final S3 transition, P(s3(3)='{s3_target}' | s2(2)='{s2_at_t2}'), is {prob_s3_final_step}.")

    # --- Step 3: Combine probabilities for the final answer ---
    print("\nStep 3: Combine the probabilities.")
    print("The total probability is the product of the S1 path probability and the final S3 transition probability.")
    final_prob = prob_s1_path * prob_s3_final_step

    print(f"\nFinal Equation: P({s1_initial}->{s1_at_t1}) * P({s1_at_t1}->{s1_at_t2}) * P({s1_at_t2}->{s1_target}) * P(s3={s3_target}|s2={s2_at_t2})")
    print(f"= {p_s1_step1} * {p_s1_step2} * {p_s1_step3} * {prob_s3_final_step}")
    print(f"= {final_prob:.4f}")
    
if __name__ == "__main__":
    solve_system_probability()
<<<0.0756>>>