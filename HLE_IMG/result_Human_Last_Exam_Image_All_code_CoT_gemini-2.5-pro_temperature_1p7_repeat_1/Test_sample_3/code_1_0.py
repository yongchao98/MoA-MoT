import numpy as np

def solve_state_values():
    """
    Calculates the state values for the given Markov Reward Process
    using Value Iteration and displays the results.
    """
    # 1. Model the Problem
    states = {
        0: 'Wake Up',
        1: 'Exercise',
        2: 'Browse Social Media',
        3: 'Work',
        4: 'Watch a Movie',
        5: 'Sleep'
    }
    num_states = len(states)

    # Rewards R(s)
    R = np.array([1.0, 3.0, -2.0, 2.0, 1.0, 0.0])

    # Discount factor gamma
    gamma = 0.2

    # Transition probability matrix P(s, s')
    P = np.zeros((num_states, num_states))

    # From Wake Up (S0)
    P[0, 1] = 0.2  # to Exercise
    P[0, 2] = 0.5  # to Browse Social Media
    P[0, 3] = 0.1  # to Work
    P[0, 4] = 0.1  # to Watch a Movie
    P[0, 5] = 0.1  # to Sleep

    # From Exercise (S1)
    P[1, 2] = 0.2  # to Browse Social Media
    P[1, 3] = 0.5  # to Work
    P[1, 5] = 0.3  # to Sleep

    # From Browse Social Media (S2)
    P[2, 1] = 0.4  # to Exercise
    P[2, 3] = 0.6  # to Work

    # From Work (S3) - Probabilities sum to 1.1, so they are normalized.
    work_transitions = np.array([0.1, 0.1, 0.3, 0.2, 0.4])
    work_transitions_normalized = work_transitions / np.sum(work_transitions)
    P[3, 0] = work_transitions_normalized[0]  # to Wake Up
    P[3, 1] = work_transitions_normalized[1]  # to Exercise
    P[3, 2] = work_transitions_normalized[2]  # to Browse Social Media
    P[3, 4] = work_transitions_normalized[3]  # to Watch a Movie
    P[3, 5] = work_transitions_normalized[4]  # to Sleep

    # From Watch a Movie (S4)
    P[4, 3] = 0.1  # to Work
    P[4, 4] = 0.7  # to Watch a Movie (self-loop)
    P[4, 5] = 0.2  # to Sleep

    # State 5 (Sleep) is terminal. Its outgoing probabilities are all 0.

    # 2. Value Iteration Algorithm
    V = np.zeros(num_states)
    tolerance = 1e-12

    while True:
        V_old = V.copy()
        for s in range(num_states):
            # The value of the terminal state 'Sleep' is its reward, which is 0.
            if s == 5:
                V[s] = R[s]
                continue
            
            # Bellman update using dot product for the summation
            V[s] = R[s] + gamma * np.dot(P[s, :], V_old)

        # Check for convergence
        if np.max(np.abs(V - V_old)) < tolerance:
            break

    # 3. Present the Results
    print("--- Final State Values ---")
    for s in range(num_states):
        print(f"V({states[s]:<20}) = {V[s]:.8f}")
    
    print("\n--- Bellman Equations with Final Values ---")
    for s in range(num_states):
        if s == 5: # Special case for the terminal state
            print(f"For state '{states[s]}':")
            print(f"{V[s]:.8f} = {R[s]}")
            print("-" * 20)
            continue
            
        # Build the string for the full equation with final numbers
        sum_terms = []
        for s_prime in range(num_states):
            if P[s, s_prime] > 0:
                # Store the string "P(s'|s) * V(s')"
                sum_terms.append(f"{P[s, s_prime]:.9f} * {V[s_prime]:.8f}")
        
        # Assemble the full equation string
        equation_str = f"{V[s]:.8f} = {R[s]} + {gamma} * ({' + '.join(sum_terms)})"
        
        print(f"For state '{states[s]}':")
        print(equation_str)
        print("-" * 20)

    # 4. Final Answer
    most_promising_idx = np.argmax(V)
    most_promising_state_value = V[most_promising_idx]

    print(f"\nThe most promising state is '{states[most_promising_idx]}' with a value of {most_promising_state_value:.8f}.")
    print(f"<<<{most_promising_state_value:.8f}>>>")

if __name__ == '__main__':
    solve_state_values()