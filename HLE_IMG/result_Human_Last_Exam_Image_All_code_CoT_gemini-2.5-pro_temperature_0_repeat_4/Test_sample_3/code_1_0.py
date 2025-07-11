import numpy as np

def solve_working_life_problem():
    """
    This function solves the working life problem using the Bellman equation.
    It calculates the state values for a given Markov Reward Process and
    identifies the value of the most promising state.
    """
    # 1. Define the states, rewards, and discount factor
    states = ["Wake Up", "Exercise", "Browse Social Media", "Work", "Watch a Movie", "Sleep"]
    # Rewards R(s) for each state
    R = np.array([1, 3, -2, 2, 1, 0], dtype=float)
    # Discount factor gamma
    gamma = 0.2

    # 2. Define the transition probability matrix P[i, j] = P(s_j | s_i)
    P = np.array([
        #  WU   Ex   BSM   Work  Movie Sleep
        [0.0, 0.2, 0.5, 0.1, 0.1, 0.1],  # From Wake Up
        [0.0, 0.0, 0.2, 0.5, 0.0, 0.3],  # From Exercise
        [0.0, 0.4, 0.0, 0.6, 0.0, 0.0],  # From Browse Social Media
        [0.1, 0.1, 0.3, 0.0, 0.2, 0.3],  # From Work
        [0.0, 0.0, 0.0, 0.1, 0.7, 0.2],  # From Watch a Movie
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]   # From Sleep (Terminal state)
    ], dtype=float)

    num_states = len(states)

    # 3. Solve the Bellman equation system: (I - gamma * P) * V = R
    I = np.identity(num_states)
    A = I - gamma * P
    V = np.linalg.solve(A, R)

    # 4. Print the final equations with the calculated values
    print("State-Value Bellman Equations (V(s) = R(s) + gamma * sum[P(s'|s) * V(s')])\n")
    for i in range(num_states):
        # Print the equation with numbers
        value_str = f"V({states[i]:<19}) = {R[i]:>4.1f} + {gamma} * ("
        value_terms = []
        for j in range(num_states):
            if P[i, j] > 0:
                value_terms.append(f"{P[i, j]} * {V[j]:.8f}")
        
        if not value_terms:
            value_str += "0)"
        else:
            # Format for better alignment
            value_str += " + ".join(value_terms) + ")"
        
        print(f"Equation for state '{states[i]}':")
        print(value_str)
        
        # Print the result of the equation
        check_val = R[i] + gamma * np.dot(P[i], V)
        print(f"{V[i]:.8f} = {check_val:.8f}\n")

    # 5. Find and print the state value of the most promising state
    max_value = np.max(V)
    most_promising_state_index = np.argmax(V)
    most_promising_state = states[most_promising_state_index]

    print("-" * 70)
    print(f"The most promising state is '{most_promising_state}' with a value of {max_value:.8f}.")
    print("\nState value of the most promising state:")
    print(f"{max_value:.8f}")

if __name__ == "__main__":
    solve_working_life_problem()