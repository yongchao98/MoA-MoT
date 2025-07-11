def solve_system_probability():
    """
    Calculates the probability of the system reaching state (D, Z, N) after 3 transitions
    from the initial state (A, X, M).
    """

    # Define the transition probabilities for S1
    # P(S1_t+1 | S1_t)
    p_s1 = {
        'A': {'B': 0.3, 'C': 0.5, 'D': 0.2},
        'B': {'A': 0.4, 'C': 0.6},
        'C': {'A': 0.1, 'B': 0.2, 'D': 0.7},
        'D': {'A': 0.8, 'B': 0.2}
    }

    # Define the S2 -> S3 transition probabilities
    # P(S3_t+1 | S2_t)
    p_s3_given_s2 = {
        'X': {'M': 0.7, 'N': 0.3},
        'Y': {'M': 0.4, 'N': 0.6},
        'Z': {'M': 0.2, 'N': 0.8}
    }

    # As determined by the step-by-step analysis, there is only one valid path for S1:
    # t=0: A -> t=1: B -> t=2: C -> t=3: D
    # The corresponding S2 state at t=2 is Y, which dictates the final S3 transition.

    # 1. Probability of S1 transitioning from A to B at step 1
    p_A_to_B = p_s1['A']['B']

    # 2. Probability of S1 transitioning from B to C at step 2
    p_B_to_C = p_s1['B']['C']

    # 3. Probability of S1 transitioning from C to D at step 3
    p_C_to_D = p_s1['C']['D']

    # 4. The S2 state at t=2 is determined by S1 at t=1 (which is B).
    # The mapping is B -> Y. The final S3 transition probability depends on this state.
    # Probability of S3 transitioning to N given S2 is Y.
    p_N_given_Y = p_s3_given_s2['Y']['N']

    # The total probability is the product of these individual probabilities.
    total_probability = p_A_to_B * p_B_to_C * p_C_to_D * p_N_given_Y

    print("The final probability is the product of the probabilities of each necessary step in the unique valid path.")
    print("The S1 state sequence must be A -> B -> C -> D for the final state to be (D, Z, N).")
    print(f"P(S1: A->B) = {p_A_to_B}")
    print(f"P(S1: B->C) = {p_B_to_C}")
    print(f"P(S1: C->D) = {p_C_to_D}")
    print("The S3 transition to N at the final step depends on the S2 state at t=2, which is Y (determined by S1 being B at t=1).")
    print(f"P(S3->N | S2=Y) = {p_N_given_Y}")
    print("\nFinal Calculation:")
    print(f"{p_A_to_B} * {p_B_to_C} * {p_C_to_D} * {p_N_given_Y} = {total_probability}")


solve_system_probability()
<<<0.0756>>>