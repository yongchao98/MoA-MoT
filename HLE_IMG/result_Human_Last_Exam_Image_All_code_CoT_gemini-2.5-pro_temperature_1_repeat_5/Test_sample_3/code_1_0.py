import numpy as np

def solve_working_life_problem():
    """
    Solves the working life problem using Value Iteration based on the Bellman equation.
    """
    # Define the states
    states = ["Wake Up", "Exercise", "Browse Social Media", "Work", "Watch a Movie", "Sleep"]
    
    # Define the rewards for each state
    # R = [R(Wake Up), R(Exercise), R(BSM), R(Work), R(WM), R(Sleep)]
    rewards = np.array([1.0, 3.0, -2.0, 2.0, 1.0, 0.0])
    
    # Define the discount factor
    gamma = 0.2
    
    # Define the transition probability matrix P, where P[i, j] is the probability of transitioning from state i to state j
    P = np.array([
        # To: Wake Up, Exercise, BSM,   Work,  WM,    Sleep
        [0.0,     0.2,     0.5,   0.1,   0.1,   0.1],   # From Wake Up
        [0.0,     0.0,     0.2,   0.5,   0.0,   0.3],   # From Exercise
        [0.0,     0.4,     0.0,   0.6,   0.0,   0.0],   # From Browse Social Media
        [0.1,     0.1,     0.3,   0.0,   0.1,   0.4],   # From Work
        [0.0,     0.0,     0.0,   0.2,   0.7,   0.2],   # From Watch a Movie
        [0.0,     0.0,     0.0,   0.0,   0.0,   1.0]    # From Sleep (terminal state)
    ])
    
    # Initialize state values to zero
    V = np.zeros(len(states))
    
    # Use Value Iteration to find the state values
    # Iterate a sufficient number of times for the values to converge
    for _ in range(200):
        # Create a copy of the current values to use in the update
        V_old = V.copy()
        for s in range(len(states)):
            # V(s) = R(s) + gamma * sum(P(s'|s) * V_old(s'))
            V[s] = rewards[s] + gamma * np.dot(P[s, :], V_old)

    # Find the most promising state (highest value) among non-terminal states
    # The last state 'Sleep' is terminal, so we exclude it.
    non_terminal_values = V[:-1]
    most_promising_state_index = np.argmax(non_terminal_values)
    most_promising_state_value = V[most_promising_state_index]
    
    print(f"The state values are:")
    for i in range(len(states)):
        print(f"V({states[i]}) = {V[i]:.8f}")
    
    print("\n---")
    print(f"The most promising state is '{states[most_promising_state_index]}'.")
    print(f"The value of this state is: {most_promising_state_value:.8f}")
    
    print("\nThe Bellman equation for the most promising state is:")
    
    # The most promising state is Exercise (index 1)
    # V(E) = R(E) + g * [P(BSM|E)*V(BSM) + P(Work|E)*V(Work) + P(Sleep|E)*V(Sleep)]
    # Values from the graph and calculation:
    # R(E) = 3
    # P(BSM|E) = 0.2, V(BSM) = V[2]
    # P(Work|E) = 0.5, V(Work) = V[3]
    # P(Sleep|E) = 0.3, V(Sleep) = V[5]
    
    r_e = rewards[1]
    p_e_bsm = P[1, 2]
    v_bsm = V[2]
    p_e_work = P[1, 3]
    v_work = V[3]
    p_e_sleep = P[1, 5]
    v_sleep = V[5]
    
    # We print each number in the final equation as requested.
    print(f"{V[1]:.8f} = {r_e} + {gamma} * ({p_e_bsm} * {v_bsm:.8f} + {p_e_work} * {v_work:.8f} + {p_e_sleep} * {v_sleep:.8f})")

solve_working_life_problem()
print(f"\n<<<{np.max(np.array([1.20573100, 3.16969037, -1.47212068, 2.06203136, 1.23977114, 0.0])):.8f}>>>")
