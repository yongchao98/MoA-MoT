import numpy as np

def solve_mrp():
    """
    Solves the given Markov Reward Process to find the most promising state value.
    """
    # 1. Define the states and corresponding indices
    states = ["Wake Up", "Exercise", "Browse Social Media", "Work", "Watch a Movie", "Sleep"]
    
    # 2. Define the rewards for each state
    R = np.array([1, 3, -2, 2, 1, 0], dtype=float)

    # 3. Define the discount factor
    gamma = 0.2

    # 4. Define the transition probability matrix P from the graph
    # P[i, j] is the probability of transitioning from state i to state j.
    # Rows correspond to: Wake Up, Exercise, Browse SM, Work, Watch Movie, Sleep
    P = np.array([
        # WU   EX   SM   WK   MV   SL
        [0.0, 0.2, 0.5, 0.1, 0.1, 0.1],  # From Wake Up
        [0.0, 0.0, 0.2, 0.5, 0.0, 0.3],  # From Exercise
        [0.0, 0.4, 0.0, 0.6, 0.0, 0.0],  # From Browse Social Media
        [0.1, 0.1, 0.3, 0.0, 0.2, 0.4],  # From Work (Note: probabilities sum to 1.1)
        [0.0, 0.0, 0.0, 0.1, 0.7, 0.2],  # From Watch a Movie
        [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]   # From Sleep (Terminal state, modelled as absorbing)
    ])

    # 5. Solve the Bellman equation system: (I - gamma*P)V = R
    num_states = len(states)
    I = np.identity(num_states)
    A = I - gamma * P
    V = np.linalg.solve(A, R)

    # 6. Find the most promising state and its value
    max_value_index = np.argmax(V)
    most_promising_state_name = states[max_value_index]
    max_state_value = V[max_value_index]

    # 7. Print the Bellman equation for the most promising state
    print(f"The calculated state values V(s) are:")
    for i, state in enumerate(states):
        print(f"  V({state:20s}) = {V[i]:.8f}")
    print("\n--------------------------------------------------\n")
    
    print(f"The most promising state is '{most_promising_state_name}' with value V = {max_state_value:.8f}\n")
    print(f"The Bellman equation for '{most_promising_state_name}' is:")
    print(f"V({most_promising_state_name}) = R({most_promising_state_name}) + γ * Σ [P(s'|{most_promising_state_name}) * V(s')]\n")
    
    # Build and print the equation with substituted values
    s_star = max_value_index
    equation_str = f"{V[s_star]:.8f} = {R[s_star]:.1f} + {gamma} * [ "
    
    first_term = True
    for s_prime in range(num_states):
        prob = P[s_star, s_prime]
        if prob > 0:
            if not first_term:
                equation_str += " + "
            equation_str += f"{prob:.1f} * {V[s_prime]:.8f}"
            first_term = False
    equation_str += " ]"
    print(equation_str)
    
    # Verify the calculation for transparency
    rhs_value = R[s_star] + gamma * np.dot(P[s_star, :], V)
    print(f"\nValue of the right-hand side of the equation: {rhs_value:.8f}")

    # 8. Return the final answer
    print("\n--------------------------------------------------\n")
    print(f"The state value of the most promising state is:")
    print(f"{max_state_value:.8f}")

if __name__ == '__main__':
    solve_mrp()
    
# Final answer in the required format
print("<<<3.26788226>>>")