import numpy as np

def solve_bellman_equation():
    """
    Solves the state-value Bellman equation for the given Markov Decision Process.
    """
    # 1. Define the states, rewards, and discount factor
    states = ["Wake Up", "Exercise", "Browse Social Media", "Work", "Watch a Movie", "Sleep"]
    n_states = len(states)
    rewards = np.array([1.0, 3.0, -2.0, 2.0, 1.0, 0.0])
    gamma = 0.2

    # 2. Define the transition probability matrix (P)
    # The probabilities from the "Work" state (index 3) in the graph sum to 1.1.
    # We normalize them to ensure they form a valid probability distribution.
    p_work_raw = np.array([0.1, 0.1, 0.3, 0.0, 0.2, 0.4])
    p_work_normalized = p_work_raw / np.sum(p_work_raw)

    P = np.array([
        # Transitions from Wake Up (0)
        [0.0, 0.2, 0.5, 0.1, 0.1, 0.1],
        # Transitions from Exercise (1)
        [0.0, 0.0, 0.2, 0.5, 0.0, 0.3],
        # Transitions from Browse Social Media (2)
        [0.0, 0.4, 0.0, 0.6, 0.0, 0.0],
        # Transitions from Work (3) - Normalized
        p_work_normalized,
        # Transitions from Watch a Movie (4)
        [0.0, 0.0, 0.0, 0.1, 0.7, 0.2],
        # Transitions from Sleep (5) - Terminal state
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    ])

    # 3. Solve the Bellman equation system: (I - gamma * P) * V = R
    I = np.identity(n_states)
    A = I - gamma * P
    V = np.linalg.solve(A, rewards)

    # 4. Print the calculated state values
    print("Calculated State Values (V):")
    for i, state in enumerate(states):
        print(f"V({state:.<20}): {V[i]:.8f}")
    print("-" * 40)

    # 5. Identify the most promising state and its value
    max_idx = np.argmax(V)
    max_value = V[max_idx]
    max_state_name = states[max_idx]

    print(f"The most promising state is '{max_state_name}' with a value of {max_value:.8f}.\n")

    # 6. Print the Bellman equation for the most promising state with all numbers
    print(f"The Bellman equation for the most promising state '{max_state_name}' is:")
    print(f"V({max_state_name}) = R({max_state_name}) + γ * Σ [ P(s'|{max_state_name}) * V(s') ]\n")
    
    equation_str = f"{V[max_idx]:.8f} = {rewards[max_idx]} + {gamma} * [ "
    
    # Build the summation part of the equation string
    summation_terms = []
    for j, state_prime in enumerate(states):
        prob = P[max_idx, j]
        if prob > 0:
            summation_terms.append(f"{prob:.8f} * {V[j]:.8f}")
    
    equation_str += " + ".join(summation_terms)
    equation_str += " ]"
    
    print("Substituting the calculated values:")
    print(equation_str)

    # Return the final answer in the required format
    return f"<<<{max_value:.8f}>>>"

# Execute the function and print the final answer
final_answer = solve_bellman_equation()
print(final_answer)
