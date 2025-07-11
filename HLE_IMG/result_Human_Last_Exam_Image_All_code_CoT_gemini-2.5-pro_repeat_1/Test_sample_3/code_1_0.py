import numpy as np

def solve_state_values():
    """
    Calculates the state values for the given Markov Reward Process
    using the Bellman equation.
    """
    # 1. Define the states, rewards, and discount factor
    states = ["Wake Up", "Exercise", "Browse Social Media", "Work", "Watch a Movie", "Sleep"]
    # Rewards for non-terminal states + terminal state
    R = np.array([1.0, 3.0, -2.0, 2.0, 1.0, 0.0])
    gamma = 0.2

    # 2. Define the transition probability matrix P
    # Rows correspond to 'from state', columns to 'to state'
    # Order: Wake Up, Exercise, BSM, Work, Movie, Sleep
    P = np.array([
        #   WU,   Ex,  BSM, Work, Movie, Sleep
        [ 0.0,  0.2,  0.5,  0.1,   0.1,   0.1 ],  # From Wake Up
        [ 0.0,  0.0,  0.2,  0.5,   0.0,   0.3 ],  # From Exercise
        [ 0.0,  0.4,  0.0,  0.6,   0.0,   0.0 ],  # From Browse Social Media
        [ 0.1,  0.1,  0.3,  0.0,   0.2,   0.3 ],  # From Work
        [ 0.0,  0.0,  0.0,  0.1,   0.7,   0.2 ],  # From Watch a Movie
        [ 0.0,  0.0,  0.0,  0.0,   0.0,   1.0 ]   # From Sleep (terminal)
    ])

    # 3. Separate terminal and non-terminal states
    num_states = len(states)
    non_terminal_indices = [0, 1, 2, 3, 4]
    num_non_terminal = len(non_terminal_indices)

    P_nt = P[np.ix_(non_terminal_indices, non_terminal_indices)]
    R_nt = R[non_terminal_indices]
    
    # The rewards vector also needs to account for transitions from non-terminal
    # states to the terminal one. The value of the terminal state V(Sleep) is R(Sleep) = 0.
    # The term gamma * P(s'|s) * V(s') for s'=Sleep is therefore 0.
    # So we only need R_nt for the right-hand side of the equation.

    # 4. Solve the system of linear equations: (I - gamma * P_nt) * V = R_nt
    I = np.identity(num_non_terminal)
    A = I - gamma * P_nt
    
    # Solve for the values of the non-terminal states
    V_nt = np.linalg.solve(A, R_nt)

    # Combine with the value of the terminal state
    V = np.zeros(num_states)
    V[non_terminal_indices] = V_nt
    V[5] = R[5] # V(Sleep) = 0

    print("--- Calculated State Values ---")
    for i in range(num_states):
        print(f"V({states[i]}) = {V[i]:.8f}")
    print("-" * 30)

    # 5. Find the most promising state
    max_value_index = np.argmax(V)
    most_promising_state_name = states[max_value_index]
    most_promising_state_value = V[max_value_index]

    print(f"The most promising state is '{most_promising_state_name}' with a value of {most_promising_state_value:.8f}\n")

    # 6. Display the Bellman equation for the most promising state
    print(f"--- Bellman Equation for '{most_promising_state_name}' ---")
    
    # Build the equation string
    equation_str = f"V({most_promising_state_name}) = R({most_promising_state_name}) + γ * Σ[P(s'|{most_promising_state_name}) * V(s')]"
    
    # Substitute values
    reward = R[max_value_index]
    summation_terms = []
    summation_value = 0
    
    for j in range(num_states):
        prob = P[max_value_index, j]
        if prob > 0:
            summation_terms.append(f"{prob} * {V[j]:.8f}")
            summation_value += prob * V[j]

    print(f"{most_promising_state_value:.8f} = {reward} + {gamma} * ({' + '.join(summation_terms)})")
    print(f"{most_promising_state_value:.8f} = {reward} + {gamma} * ({summation_value:.8f})")
    print(f"{most_promising_state_value:.8f} = {reward + gamma * summation_value:.8f}")

    # Return the final answer as a formatted string
    return f"{most_promising_state_value:.8f}"

if __name__ == '__main__':
    final_answer = solve_state_values()
    print(f"\n<<<_START_END_>>>\n{final_answer}\n<<<!!!>>>")
