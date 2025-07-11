import numpy as np

def solve_working_life_problem():
    """
    Solves the state-value Bellman equation for the given Markov Reward Process.
    """
    # 1. Define the model parameters
    states = ['Wake Up', 'Exercise', 'Browse Social Media', 'Work', 'Watch a Movie', 'Sleep']
    # Rewards for each state
    rewards = np.array([1, 3, -2, 2, 1, 0])
    # Discount factor
    gamma = 0.2

    # Transition probability matrix P[i, j] = P(s_j | s_i)
    # The order of states is as defined in the 'states' list.
    P = np.array([
        #           WU,  EX,  BSM, WORK, MOVIE, SLEEP
        [0.0, 0.2, 0.5, 0.1, 0.1, 0.1],   # Transitions from Wake Up
        [0.0, 0.0, 0.2, 0.5, 0.0, 0.3],   # Transitions from Exercise
        [0.0, 0.4, 0.0, 0.6, 0.0, 0.0],   # Transitions from Browse Social Media
        [0.1, 0.1, 0.3, 0.0, 0.2, 0.4],   # Transitions from Work
        [0.0, 0.0, 0.0, 0.1, 0.7, 0.2],   # Transitions from Watch a Movie
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]    # Transitions from Sleep (terminal state)
    ])

    # 2. Separate non-terminal and terminal states for solving
    # Non-terminal states are the first 5 (indices 0 to 4)
    # Terminal state 'Sleep' is index 5
    num_non_terminal_states = 5
    P_nt = P[0:num_non_terminal_states, 0:num_non_terminal_states]
    R_nt = rewards[0:num_non_terminal_states]

    # 3. Solve the system of linear equations: (I - gamma * P_nt) * V = R_nt
    I = np.identity(num_non_terminal_states)
    A = I - gamma * P_nt
    
    # Calculate the values for the non-terminal states
    V_nt = np.linalg.solve(A, R_nt)

    # Combine with terminal state value to get the full value vector
    V = np.append(V_nt, rewards[5])

    # 4. Find the most promising state and its value
    max_value_index = np.argmax(V)
    most_promising_state_name = states[max_value_index]
    most_promising_state_value = V[max_value_index]

    # 5. Output the results, including the detailed Bellman equation for the best state
    print("Calculated State Values:")
    for i in range(len(states)):
        print(f"V({states[i]}) = {V[i]:.8f}")
    print("\n--------------------------------------------------")
    print(f"The most promising state is '{most_promising_state_name}'.")
    
    # Construct and print the Bellman equation with all numbers
    print("\nThe Bellman equation for this state with all numerical values is:")
    
    reward = rewards[max_value_index]
    equation_rhs_terms = []
    
    for j in range(len(states)):
        prob = P[max_value_index, j]
        if prob > 0:
            term = f"{prob} * {V[j]:.8f}"
            equation_rhs_terms.append(term)
            
    equation_str = f"{most_promising_state_value:.8f} = {reward} + {gamma} * ({' + '.join(equation_rhs_terms)})"
    print(equation_str)
    
    print("\n--------------------------------------------------")
    print("The state value of the most promising state is:")
    print(f"{most_promising_state_value:.8f}")
    
    # Final answer in the required format
    print(f"\n<<<{most_promising_state_value:.8f}>>>")

if __name__ == '__main__':
    solve_working_life_problem()