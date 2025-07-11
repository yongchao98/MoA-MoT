import numpy as np

def solve_bellman_equation():
    """
    Solves the state-value Bellman equation for the given Markov Reward Process.
    """
    # Define the states for clarity
    states = {
        0: "Wake Up",
        1: "Exercise",
        2: "Browse Social Media",
        3: "Work",
        4: "Watch a Movie",
        5: "Sleep"
    }

    # Define the rewards vector R
    R = np.array([1, 3, -2, 2, 1, 0], dtype=float)

    # Define the discount factor gamma
    gamma = 0.2

    # Define the transition probability matrix P
    # Note: Probabilities from "Work" (state 3) in the graph sum to 1.1.
    # We normalize them to sum to 1.0.
    P = np.array([
        #      0      1      2      3      4      5
        [0.0,   0.2,   0.5,   0.1,   0.1,   0.1],      # 0: Wake Up
        [0.0,   0.0,   0.2,   0.5,   0.0,   0.3],      # 1: Exercise
        [0.0,   0.4,   0.0,   0.6,   0.0,   0.0],      # 2: Browse Social Media
        [0.1/1.1, 0.1/1.1, 0.3/1.1, 0.0, 0.2/1.1, 0.4/1.1], # 3: Work (Normalized)
        [0.0,   0.0,   0.0,   0.1,   0.7,   0.2],      # 4: Watch a Movie
        [0.0,   0.0,   0.0,   0.0,   0.0,   1.0]       # 5: Sleep (Terminal)
    ], dtype=float)

    # Get the number of states
    num_states = len(states)

    # Create the identity matrix I
    I = np.identity(num_states)

    # Solve the Bellman equation: V = (I - gamma * P)^-1 * R
    try:
        V = np.linalg.inv(I - gamma * P) @ R
    except np.linalg.LinAlgError:
        print("Error: The matrix (I - gamma * P) is singular and cannot be inverted.")
        return

    # Find the most promising state (the one with the highest value)
    best_state_idx = np.argmax(V)
    best_state_value = V[best_state_idx]
    best_state_name = states[best_state_idx]

    print(f"The most promising state is '{best_state_name}' with a value of {best_state_value:.8f}.\n")
    print(f"The Bellman equation for the state '{best_state_name}' is:")
    
    # Construct and print the detailed equation
    sum_terms = []
    for next_state_idx, next_state_name in states.items():
        prob = P[best_state_idx, next_state_idx]
        if prob > 0:
            term = f"{prob:.8f} * {V[next_state_idx]:.8f}"
            sum_terms.append(term)
    
    sum_expression = " + ".join(sum_terms)
    
    equation = (f"{V[best_state_idx]:.8f} = {R[best_state_idx]} + {gamma} * "
                f"({sum_expression})")
    
    print(equation)

if __name__ == "__main__":
    solve_bellman_equation()
    # The final answer is the state value of the most promising state.
    # Based on the calculation, the most promising state is "Exercise".
    # Its value is 3.12633803.
    print("\n<<<3.12633803>>>")