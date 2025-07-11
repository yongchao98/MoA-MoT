import numpy as np

def solve_state_values():
    """
    Solves the Bellman equation for the given working life problem.
    """
    # --- 1. Define the components of the Markov Reward Process ---
    
    # States mapping for clarity
    state_names = ["Wake Up", "Exercise", "Browse Social Media", "Work", "Watch a Movie", "Sleep"]
    num_states = len(state_names)
    
    # Rewards vector R
    rewards = np.array([1, 3, -2, 2, 1, 0], dtype=float)
    
    # Discount factor gamma
    gamma = 0.2
    
    # Transition probability matrix P
    P = np.zeros((num_states, num_states), dtype=float)
    
    # Transitions from Wake Up (State 0)
    P[0] = [0.0, 0.2, 0.5, 0.1, 0.1, 0.1]
    
    # Transitions from Exercise (State 1)
    P[1] = [0.0, 0.0, 0.2, 0.5, 0.0, 0.3]
    
    # Transitions from Browse Social Media (State 2)
    P[2] = [0.0, 0.4, 0.0, 0.6, 0.0, 0.0]
    
    # Transitions from Work (State 3)
    # The sum of probabilities from Work in the graph is 1.1. We normalize them.
    work_transitions = np.array([0.1, 0.1, 0.3, 0.0, 0.2, 0.4])
    P[3] = work_transitions / np.sum(work_transitions)
    
    # Transitions from Watch a Movie (State 4)
    P[4] = [0.0, 0.0, 0.0, 0.1, 0.7, 0.2]
    
    # Transitions from Sleep (State 5) - Terminal state
    P[5] = [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]

    # --- 2. Solve the Bellman Equation ---
    
    # We solve the system of linear equations (I - gamma*P)V = R
    I = np.identity(num_states)
    A = I - gamma * P
    
    # Solve for V
    V = np.linalg.solve(A, rewards)

    # --- 3. Display the results ---
    
    print("--- State-Value Bellman Equation Solution ---\n")
    print("The calculated state-values V(s) are:")
    for i in range(num_states):
        print(f"  V({state_names[i]:<20}) = {V[i]:.8f}")
    
    print("\n--- Verification of Bellman Equations ---")
    print("V(s) = R(s) + γ * Σ [P(s'|s) * V(s')]\n")

    for i in range(num_states):
        print(f"For state '{state_names[i]}':")
        # Calculate the right-hand side of the equation to verify
        rhs = rewards[i] + gamma * np.sum(P[i] * V)
        
        # Build the equation string with all numbers
        equation_str = f"  {V[i]:.8f} = {rewards[i]} + {gamma} * ("
        first_term = True
        for j in range(num_states):
            if P[i, j] > 0:
                if not first_term:
                    equation_str += " + "
                # We show the full precision probability used in the calculation
                equation_str += f"{P[i, j]} * {V[j]:.8f}"
                first_term = False
        equation_str += f") = {rhs:.8f}"
        print(equation_str)
        print("-" * 20)
        
    # --- 4. Find the most promising state ---
    
    most_promising_state_index = np.argmax(V)
    highest_value = V[most_promising_state_index]
    most_promising_state_name = state_names[most_promising_state_index]
    
    print(f"\nThe most promising state is '{most_promising_state_name}'.")
    print(f"The state value of this state is {highest_value:.8f}.")

    return highest_value

if __name__ == '__main__':
    result = solve_state_values()
    # The final answer in the requested format
    print(f"\n<<<{result:.8f}>>>")
