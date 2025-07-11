import numpy as np

def solve_mrp():
    """
    Solves the working life problem by calculating state values using the Bellman equation.
    """
    # 1. Define the model components
    states = {
        "Wake Up": 0, "Exercise": 1, "Browse Social Media": 2,
        "Work": 3, "Watch a Movie": 4, "Sleep": 5
    }
    state_names = list(states.keys())
    
    R_full = np.array([1, 3, -2, 2, 1, 0], dtype=float)
    gamma = 0.2

    # 2. Handle inconsistencies in the transition probability matrix P
    # Probabilities from "Watch a Movie" are normalized (0.2 + 0.2 + 0.7 = 1.1)
    p_watch_work = 0.2
    p_watch_sleep = 0.2
    p_watch_self = 0.7
    total_p_watch = p_watch_work + p_watch_sleep + p_watch_self

    P_full = np.array([
        # Transitions from ... to -> WakeUp, Exercise, Browse, Work, Movie, Sleep
        [0.0, 0.2, 0.5, 0.1, 0.1, 0.1],            # Wake Up
        [0.0, 0.0, 0.2, 0.5, 0.0, 0.3],            # Exercise
        [0.0, 0.4, 0.0, 0.6, 0.0, 0.0],            # Browse Social Media
        [0.1, 0.1, 0.3, 0.0, 0.1, 0.4],            # Work
        [0.0, 0.0, 0.0, p_watch_work / total_p_watch, p_watch_self / total_p_watch, p_watch_sleep / total_p_watch], # Watch a Movie
        [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]             # Sleep (Terminal state)
    ], dtype=float)

    # 3. Separate non-terminal states for solving the system of equations
    # "Sleep" is a terminal state, its value V(Sleep) is fixed at its reward, 0.
    num_non_terminal = 5
    P_nt = P_full[:num_non_terminal, :num_non_terminal]
    R_nt = R_full[:num_non_terminal]
    
    # 4. Set up and 5. solve the Bellman equation system: (I - γ * P_nt) * V_nt = R_nt
    I = np.identity(num_non_terminal)
    A = I - gamma * P_nt
    V_nt = np.linalg.solve(A, R_nt)

    # Combine results with the terminal state's value
    V_full = np.append(V_nt, R_full[states["Sleep"]])

    # 6. Find the most promising state (highest V(s))
    max_val_idx = np.argmax(V_full)
    max_val = V_full[max_val_idx]
    max_state_name = state_names[max_val_idx]
    
    # 7. Output the results and the final equation
    print("The state values V(s) are calculated by solving the Bellman equation: V = R + γPV")
    print("\nCalculated state values for all states:")
    for i in range(len(V_full)):
        print(f"V({state_names[i]:<20}) = {V_full[i]:.8f}")
    
    print(f"\nThe most promising state is '{max_state_name}' with a value of {max_val:.8f}.")

    print("\nThe Bellman equation for this state is:")
    print(f"V({max_state_name}) = R({max_state_name}) + γ * Σ [ P(s' | {max_state_name}) * V(s') ]")
    
    # Construct the detailed equation string
    equation_str = f"{max_val:.8f} = {R_full[max_val_idx]} + {gamma} * ("
    sum_terms = []
    for s_prime_idx in range(len(state_names)):
        prob = P_full[max_val_idx, s_prime_idx]
        if prob > 0:
            sum_terms.append(f"{prob:.4f} * {V_full[s_prime_idx]:.8f}")
            
    equation_str += " + ".join(sum_terms)
    equation_str += ")"
    print(equation_str)

if __name__ == '__main__':
    solve_mrp()