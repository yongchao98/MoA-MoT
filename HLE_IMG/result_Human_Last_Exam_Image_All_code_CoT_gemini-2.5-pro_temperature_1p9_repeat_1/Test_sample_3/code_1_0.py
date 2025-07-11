import numpy as np

def solve_state_values():
    """
    Calculates the state values for the given Markov Reward Process using the Bellman equation.
    """
    # Define states, rewards, and discount factor
    states = ["Wake Up", "Exercise", "Browse Social Media", "Work", "Watch a Movie"]
    # Rewards for non-terminal states [Wake Up, Exercise, BSM, Work, WaM]
    R = np.array([1, 3, -2, 2, 1], dtype=float)
    gamma = 0.2

    # As noted in the plan, the probabilities from "Work" sum to 1.1 in the graph.
    # We assume P(Sleep | Work) = 0.3 instead of 0.4, so the probabilities sum to 1.
    # P is the transition probability matrix for non-terminal states P[i, j] = P(s'=j | s=i)
    # The transitions to the terminal state 'Sleep' are implicitly handled because V(Sleep)=0.
    P = np.array([
        # To: Wake Up, Exercise, BSM,   Work,  WaM
        [0.0,     0.2,      0.5,    0.1,   0.1],   # From Wake Up
        [0.0,     0.0,      0.2,    0.5,   0.0],   # From Exercise
        [0.0,     0.4,      0.0,    0.6,   0.0],   # From Browse Social Media
        [0.1,     0.1,      0.3,    0.0,   0.2],   # From Work (with P(Work->Sleep)=0.3)
        [0.0,     0.0,      0.0,    0.1,   0.7]    # From Watch a Movie
    ])

    # The Bellman equation is V = R + gamma * P * V
    # We solve the linear system (I - gamma*P)V = R
    I = np.identity(len(states))
    A = I - gamma * P
    
    # Solve for V
    V = np.linalg.solve(A, R)
    
    print("State Values Calculation:")
    print("-------------------------")
    print(f"States (X): {tuple(states)}")
    print(f"Rewards (r_X): {list(R)}")
    print(f"Discount factor (γ): {gamma}")
    
    print("\nCalculated State Values (V) to 8 decimal points:")
    for i in range(len(states)):
        print(f"V({states[i]:<20}) = {V[i]:.8f}")
        
    # Find the most promising state (highest value)
    most_promising_idx = np.argmax(V)
    most_promising_state = states[most_promising_idx]
    max_value = V[most_promising_idx]

    print(f"\nThe most promising state is '{most_promising_state}' with a value of {max_value:.8f}.")

    print("\nThe Bellman equation for this state is:")
    print(f"V({most_promising_state}) = r({most_promising_state}) + γ * Σ [ P({most_promising_state} -> s') * V(s') ]")
    
    print("\nSubstituting the numerical values:")
    
    # Probabilities from 'Exercise'
    p_ex_bsm = 0.2
    p_ex_work = 0.5
    p_ex_sleep = 0.3  # Derived from 1 - 0.2 - 0.5
    
    r_ex = R[1]
    v_bsm = V[2]
    v_work = V[3]
    v_sleep = 0.0 # By definition
    
    print(f"{max_value:.8f} = {r_ex:.8f} + {gamma:.8f} * ({p_ex_bsm:.8f} * {v_bsm:.8f} + {p_ex_work:.8f} * {v_work:.8f} + {p_ex_sleep:.8f} * {v_sleep:.8f})")
    
    print(f"\nFinal Answer: The state value of the most promising state ('{most_promising_state}') is {max_value:.8f}")

if __name__ == "__main__":
    solve_state_values()
    # To extract the final numerical answer for the platform
    # The logic is re-run here just for the final print, avoiding global variables
    R = np.array([1, 3, -2, 2, 1], dtype=float)
    gamma = 0.2
    P = np.array([
        [0.0, 0.2, 0.5, 0.1, 0.1],
        [0.0, 0.0, 0.2, 0.5, 0.0],
        [0.0, 0.4, 0.0, 0.6, 0.0],
        [0.1, 0.1, 0.3, 0.0, 0.2],
        [0.0, 0.0, 0.0, 0.1, 0.7]
    ])
    I = np.identity(5)
    A = I - gamma * P
    V = np.linalg.solve(A, R)
    max_value = np.max(V)
    print(f"<<<{max_value:.8f}>>>")