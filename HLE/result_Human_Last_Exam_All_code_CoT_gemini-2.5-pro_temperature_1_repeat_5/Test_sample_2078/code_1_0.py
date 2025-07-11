def solve_system_probability():
    """
    Calculates the probability of the system being in state (D, Z, N) after 3 transitions.
    """

    # Define the transition probabilities for S1
    s1_transitions = {
        'A': {'B': 0.3, 'C': 0.5, 'D': 0.2},
        'B': {'A': 0.4, 'C': 0.6},
        'C': {'A': 0.1, 'B': 0.2, 'D': 0.7},
        'D': {'A': 0.8, 'B': 0.2}
    }

    # Define the S2 -> S3 probabilistic interaction
    s2_s3_transitions = {
        'X': {'M': 0.7, 'N': 0.3},
        'Y': {'M': 0.4, 'N': 0.6},
        'Z': {'M': 0.2, 'N': 0.8}
    }

    # Step 1: Identify the required S1 path probabilities.
    # To end at S1=D at t=3 from S1=C at t=2, we need P(C->D).
    # To get to S1=C at t=2 from S1=A at t=0, the only valid path is A->B->C.
    # We need P(A->B) and P(B->C).
    p_A_to_B = s1_transitions['A']['B']
    p_B_to_C = s1_transitions['B']['C']
    p_C_to_D = s1_transitions['C']['D']

    # Step 2: Identify the required S3 transition probability.
    # The S1 path A->B->C determines that S1(1)=B.
    # The state S2(2) is determined by S1(1), so S2(2) is Y.
    # The target state S3(3) is N. This transition depends on S2(2).
    # So we need the probability of S3 transitioning to N given S2 is Y.
    p_N_given_Y = s2_s3_transitions['Y']['N']

    # Step 3: Calculate the total probability.
    # The total probability is the product of the probabilities of each required event
    # in the sequence.
    # P(total) = P(A->B) * P(B->C) * P(C->D) * P(S3=N|S2(2)=Y)
    total_probability = p_A_to_B * p_B_to_C * p_C_to_D * p_N_given_Y
    
    # Print the final equation with all the numbers
    print("The probability is calculated as the product of the probabilities of the required sequence of transitions:")
    print(f"P(A->B) * P(B->C) * P(C->D) * P(S3=N | S2=Y)")
    print("The final equation is:")
    print(f"{p_A_to_B} * {p_B_to_C} * {p_C_to_D} * {p_N_given_Y} = {total_probability}")

solve_system_probability()
<<<0.0756>>>