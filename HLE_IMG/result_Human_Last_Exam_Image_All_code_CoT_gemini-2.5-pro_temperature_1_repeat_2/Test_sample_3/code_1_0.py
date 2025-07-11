import numpy as np

def solve_mrp():
    """
    Calculates the state values for the given Markov Reward Process
    using the Bellman equation.
    """
    # Define states for clarity in the output
    states = ["Wake Up", "Exercise", "Browse Social Media", "Work", "Watch a Movie", "Sleep"]
    num_states = len(states)

    # Define the rewards for each state
    # R = [r_wu, r_ex, r_bsm, r_wrk, r_wam, r_slp]
    R = np.array([1.0, 3.0, -2.0, 2.0, 1.0, 0.0])

    # Define the discount factor
    gamma = 0.2

    # Define the transition probability matrix P[i, j] = P(S' = j | S = i)
    # Note: Probabilities from "Watch a Movie" (state 4) in the diagram sum to 1.1.
    # We assume a typo and correct the self-loop from 0.7 to 0.6 (1 - 0.2 - 0.2).
    P = np.array([
        # To: WU,   EX,   BSM,  WRK,  WAM,  SLP
        [0.0,  0.2,  0.5,  0.1,  0.1,  0.1],  # From Wake Up
        [0.0,  0.0,  0.2,  0.5,  0.0,  0.3],  # From Exercise
        [0.0,  0.4,  0.0,  0.6,  0.0,  0.0],  # From Browse Social Media
        [0.1,  0.1,  0.3,  0.0,  0.1,  0.4],  # From Work
        [0.0,  0.0,  0.0,  0.2,  0.6,  0.2],  # From Watch a Movie (corrected)
        [0.0,  0.0,  0.0,  0.0,  0.0,  1.0]   # From Sleep (terminal)
    ])

    # The Bellman equation V = R + gamma * P * V can be solved as a system of linear equations:
    # (I - gamma * P) * V = R
    # We solve this for the non-terminal states, as the terminal state value is known (V_sleep = 0).

    # Indices for non-terminal states (0-4) and terminal states (5)
    nt_indices = list(range(5))
    
    # Extract sub-matrix and sub-vector for non-terminal states
    P_nt = P[np.ix_(nt_indices, nt_indices)]
    R_nt = R[nt_indices]

    # Create the identity matrix
    I = np.identity(len(nt_indices))

    # Solve the system (I - gamma * P_nt) * V_nt = R_nt
    A = I - gamma * P_nt
    V_nt = np.linalg.solve(A, R_nt)

    # Combine results into the full state-value vector V
    V = np.zeros(num_states)
    V[nt_indices] = V_nt
    # V[5] is already 0, which is the reward for the terminal state 'Sleep'

    # --- Output Results ---
    print("The solved Bellman equations for each state are:\n")

    for i in range(num_states):
        equation_lhs = f"V({states[i]})"
        
        # Handle the terminal state 'Sleep'
        if i == 5:
            print(f"{equation_lhs} = {R[i]:.8f}")
            continue

        sum_terms_str = []
        for j in range(num_states):
            prob = P[i, j]
            if prob > 0:
                sum_terms_str.append(f"{prob:.1f} * {V[j]:.8f}")
        
        equation_rhs = f"{R[i]:.8f} + {gamma} * ({' + '.join(sum_terms_str)})"
        
        print(f"{equation_lhs} = {equation_rhs} = {V[i]:.8f}\n")

    print("\n--- Summary of State Values ---")
    for i in range(num_states):
        print(f"V({states[i]}) = {V[i]:.8f}")

    most_promising_value = np.max(V)
    most_promising_state_index = np.argmax(V)

    print(f"\nThe most promising state is '{states[most_promising_state_index]}'.")
    print("\nThe state value of the most promising state is:")
    print(f"{most_promising_value:.8f}")

if __name__ == '__main__':
    solve_mrp()