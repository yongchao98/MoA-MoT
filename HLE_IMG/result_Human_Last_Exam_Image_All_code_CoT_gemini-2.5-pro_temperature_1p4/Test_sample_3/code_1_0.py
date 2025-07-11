import numpy as np

def solve_working_life_problem():
    """
    Solves the state-value Bellman equation for the given working life problem.
    """
    # 1. Model the Problem
    states = ["Wake Up", "Exercise", "Browse Social Media", "Work", "Watch a Movie", "Sleep"]
    n_states = len(states)

    # Rewards vector R
    R = np.array([1, 3, -2, 2, 1, 0], dtype=float)

    # Discount factor gamma
    gamma = 0.2

    # Transition probability matrix P
    P = np.zeros((n_states, n_states), dtype=float)

    # Transitions from Wake Up (S0)
    P[0, 1] = 0.2  # to Exercise
    P[0, 2] = 0.5  # to Browse Social Media
    P[0, 3] = 0.1  # to Work
    P[0, 4] = 0.1  # to Watch a Movie
    P[0, 5] = 0.1  # to Sleep

    # Transitions from Exercise (S1)
    P[1, 2] = 0.2  # to Browse Social Media
    P[1, 3] = 0.5  # to Work
    P[1, 5] = 0.3  # to Sleep

    # Transitions from Browse Social Media (S2)
    P[2, 1] = 0.4  # to Exercise
    P[2, 3] = 0.6  # to Work

    # Transitions from Work (S3) - Probabilities in the graph sum to 1.1, so we normalize them.
    work_probs = np.array([
        0.1,  # to Wake Up (S0)
        0.1,  # to Exercise (S1)
        0.3,  # to Browse Social Media (S2)
        0.2,  # to Watch a Movie (S4)
        0.4   # to Sleep (S5)
    ])
    work_sum = np.sum(work_probs)
    P[3, 0] = work_probs[0] / work_sum
    P[3, 1] = work_probs[1] / work_sum
    P[3, 2] = work_probs[2] / work_sum
    P[3, 4] = work_probs[3] / work_sum
    P[3, 5] = work_probs[4] / work_sum

    # Transitions from Watch a Movie (S4)
    P[4, 3] = 0.1  # to Work
    P[4, 4] = 0.7  # to itself
    P[4, 5] = 0.2  # to Sleep

    # Transitions from Sleep (S5) - Terminal state
    P[5, 5] = 1.0

    # 2 & 3. Set up and solve the Bellman equation
    # (I - gamma * P) * V = R
    I = np.identity(n_states)
    A = I - gamma * P
    V = np.linalg.solve(A, R)

    # 4. Identify the most promising state
    most_promising_state_index = np.argmax(V)
    most_promising_state_name = states[most_promising_state_index]
    max_value = V[most_promising_state_index]

    # 5. Output the results
    print("Calculated State Values (to 8 decimal points):")
    for i, state in enumerate(states):
        print(f"  V({state:20s}) = {V[i]:.8f}")

    print(f"\nThe most promising state is '{most_promising_state_name}' with a value of {max_value:.8f}.\n")
    
    print("The Bellman equation for the most promising state is:")
    print(f"V({most_promising_state_name}) = R({most_promising_state_name}) + γ * Σ [ P(s'|{most_promising_state_name}) * V(s') ]")
    
    # Constructing the equation string
    s_idx = most_promising_state_index
    equation_rhs = f"{R[s_idx]} + {gamma} * ["
    terms = []
    for s_prime_idx, prob in enumerate(P[s_idx]):
        if prob > 0:
            terms.append(f"{prob:.4f} * {V[s_prime_idx]:.8f}")
    equation_rhs += " + ".join(terms) + " ]"
    
    print(f"{V[s_idx]:.8f} = {equation_rhs}")
    
    # Print the final answer in the specified format
    print(f"\n{max_value:.8f}")

solve_working_life_problem()