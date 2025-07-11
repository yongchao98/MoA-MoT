import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

def solve_cas_probability():
    """
    Calculates the probability of the system being in state (D, Z, N) after 3 transitions.
    """
    # S1: Transition probabilities
    p_s1 = {
        'A': {'A': 0.0, 'B': 0.3, 'C': 0.5, 'D': 0.2},
        'B': {'A': 0.4, 'B': 0.0, 'C': 0.6, 'D': 0.0},
        'C': {'A': 0.1, 'B': 0.2, 'C': 0.0, 'D': 0.7},
        'D': {'A': 0.8, 'B': 0.2, 'C': 0.0, 'D': 0.0}
    }

    # S1 -> S2: Deterministic mapping S2_{t+1} = f(S1_t)
    # The mapping for A, B, C is implied. D is assumed to cycle back to X.
    # As shown in the logic, only f('B') is needed for this specific problem.
    s1_to_s2_map = {'A': 'X', 'B': 'Y', 'C': 'Z', 'D': 'X'}

    # S2 -> S3: Probabilistic transitions P(S3_{t+1} | S2_t)
    p_s3 = {
        'X': {'M': 0.7, 'N': 0.3},
        'Y': {'M': 0.4, 'N': 0.6},
        'Z': {'M': 0.2, 'N': 0.8}
    }

    initial_s1 = 'A'
    # Target S1 path structure must be: S1_0 -> S1_1 -> C -> D
    s1_2 = 'C'
    s1_3 = 'D'
    
    # Target S3 state at t=3
    s3_3 = 'N'

    total_prob = 0.0
    
    s1_states = ['A', 'B', 'C', 'D']
    
    # The final equation components
    path_components = []

    # Sum probabilities over all possible intermediate states S1_1
    for s1_1 in s1_states:
        # Probability of the S1 path segment: A -> S1_1 -> C -> D
        p_s1_path_A_to_1 = p_s1[initial_s1][s1_1]
        p_s1_path_1_to_2 = p_s1[s1_1][s1_2]
        p_s1_path_2_to_3 = p_s1[s1_2][s1_3]
        
        prob_s1_path = p_s1_path_A_to_1 * p_s1_path_1_to_2 * p_s1_path_2_to_3
        
        # If the S1 path has zero probability, it's not a valid path.
        if prob_s1_path == 0:
            continue

        # Determine the state S2_2, which influences the final S3 transition
        s2_2 = s1_to_s2_map[s1_1]

        # Probability of the S3 transition S3_2 -> N, given S2_2
        prob_s3_transition = p_s3[s2_2][s3_3]
        
        # The probability of this specific path is the product of S1 path prob and S3 transition prob
        path_prob = prob_s1_path * prob_s3_transition
        total_prob += path_prob
        
        # Store the numbers for the final equation output
        path_components.append({
            "p_A_to_s1_1": p_s1_path_A_to_1,
            "s1_1": s1_1,
            "p_s1_1_to_C": p_s1_path_1_to_2,
            "p_C_to_D": p_s1_path_2_to_3,
            "s2_2": s2_2,
            "p_N_given_s2_2": prob_s3_transition,
            "result": path_prob
        })

    print("The final state (D, Z, N) requires the S1 path to be of the form A -> S1_1 -> C -> D.")
    print("We found one valid path where S1_1 = B.")
    print("\nThe probability is the product of the S1 path probability and the conditional S3 transition probability.")
    
    # Print the equation with the numbers from the valid path
    if path_components:
        comp = path_components[0]
        print("\nFinal Equation:")
        print(f"P(total) = P(A->B) * P(B->C) * P(C->D) * P(S3=N | S2=f(B)=Y)")
        print(f"P(total) = {comp['p_A_to_s1_1']} * {comp['p_s1_1_to_C']} * {comp['p_C_to_D']} * {comp['p_N_given_s2_2']}")
        print(f"P(total) = {comp['result']}")
    else:
        print("\nNo valid paths were found. The probability is 0.")

    print(f"\nThe final probability of the system being in state (D, Z, N) after 3 transitions is: {total_prob}")
    
solve_cas_probability()

# Get the content from the buffer
output = captured_output.getvalue()
# Close the buffer
captured_output.close()
# Restore original stdout
sys.stdout = original_stdout

# Print the captured output
print(output.strip())
print(f"<<<{total_prob}>>>")