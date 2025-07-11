import numpy as np

def solve_mdp():
    """
    Solves the state-value Bellman equation for the given MDP.
    """
    # Define the states for clear output
    states = ["Wake Up", "Exercise", "Browse Social Media", "Work", "Watch a Movie", "Sleep"]
    n_states = len(states)

    # Define the rewards for each state
    # R = [R_wu, R_e, R_bsm, R_w, R_wm, R_s]
    R = np.array([1.0, 3.0, -2.0, 2.0, 1.0, 0.0])

    # Define the discount factor
    gamma = 0.2

    # Define the transition probability matrix P[s, s']
    # The probabilities from 'Work' sum to 1.1, so they are normalized.
    work_probs_sum = 0.1 + 0.1 + 0.3 + 0.2 + 0.4
    P = np.array([
        # WU    E      BSM    W      WM     S
        [0.0,   0.2,   0.5,   0.1,   0.1,   0.1],          # From Wake Up
        [0.0,   0.0,   0.2,   0.5,   0.0,   0.3],          # From Exercise
        [0.0,   0.4,   0.0,   0.6,   0.0,   0.0],          # From Browse Social Media
        [0.1/work_probs_sum, 0.1/work_probs_sum, 0.3/work_probs_sum, 0.0, 0.2/work_probs_sum, 0.4/work_probs_sum], # From Work
        [0.0,   0.0,   0.0,   0.1,   0.7,   0.2],          # From Watch a Movie
        [0.0,   0.0,   0.0,   0.0,   0.0,   0.0]           # From Sleep (Terminal State)
    ])

    # The Bellman equation is V = R + gamma * P @ V
    # We can solve this as a system of linear equations: (I - gamma * P) @ V = R
    I = np.identity(n_states)
    A = I - gamma * P
    
    # Solve for V
    V = np.linalg.solve(A, R)

    print("This problem can be modeled as a system of Bellman equations V(s) = R(s) + γ * Σ_s' P(s'|s)V(s').")
    print("Here are the equations for each state:\n")

    # Print the equations with numbers
    # V_s = ["V_wu", "V_e", "V_bsm", "V_w", "V_wm", "V_s"] # Symbolic names
    v_names = ["V(Wake Up)", "V(Exercise)", "V(Browse SM)", "V(Work)", "V(Watch Movie)", "V(Sleep)"]
    for i in range(n_states):
        equation = f"{v_names[i]} = {R[i]} + {gamma} * ("
        terms = []
        for j in range(n_states):
            if P[i, j] > 0:
                terms.append(f"{P[i, j]:.4f} * {v_names[j]}")
        if not terms:
            equation += "0)"
        else:
            equation += " + ".join(terms) + ")"
        print(equation)

    print("\nBy solving this system of linear equations, we get the following state values:")
    for i in range(n_states):
        print(f"Value of {states[i]}: {V[i]:.8f}")
        
    # Find the most promising state
    most_promising_idx = np.argmax(V)
    most_promising_value = V[most_promising_idx]
    most_promising_state = states[most_promising_idx]

    print(f"\nThe most promising state is '{most_promising_state}'.")
    print(f"The state value of the most promising state is: {most_promising_value:.8f}")

    return most_promising_value

if __name__ == '__main__':
    result = solve_mdp()
    # The final answer is wrapped in <<<>>>
    print(f"\n<<<{result:.8f}>>>")
