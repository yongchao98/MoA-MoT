import numpy as np

def solve_working_life_problem():
    """
    Calculates the state values for the working life problem using the Bellman equation.
    """
    # 1. Define states, rewards, and the discount factor
    states = ['Wake Up', 'Exercise', 'Browse Social Media', 'Work', 'Watch a Movie', 'Sleep']
    state_to_idx = {name: i for i, name in enumerate(states)}
    rewards = np.array([1, 3, -2, 2, 1, 0], dtype=float)
    gamma = 0.2

    # 2. Define the transition probability matrix P
    P = np.zeros((len(states), len(states)))

    # From 'Wake Up' (s0)
    P[state_to_idx['Wake Up'], state_to_idx['Exercise']] = 0.2
    P[state_to_idx['Wake Up'], state_to_idx['Browse Social Media']] = 0.5
    P[state_to_idx['Wake Up'], state_to_idx['Work']] = 0.1
    P[state_to_idx['Wake Up'], state_to_idx['Watch a Movie']] = 0.1
    P[state_to_idx['Wake Up'], state_to_idx['Sleep']] = 0.1

    # From 'Exercise' (s1)
    P[state_to_idx['Exercise'], state_to_idx['Browse Social Media']] = 0.2
    P[state_to_idx['Exercise'], state_to_idx['Work']] = 0.5
    P[state_to_idx['Exercise'], state_to_idx['Sleep']] = 0.3

    # From 'Browse Social Media' (s2)
    P[state_to_idx['Browse Social Media'], state_to_idx['Exercise']] = 0.4
    P[state_to_idx['Browse Social Media'], state_to_idx['Work']] = 0.6

    # From 'Work' (s3) - Probabilities normalized because they sum to 1.1 in the graph
    work_transitions = np.array([0.1, 0.1, 0.3, 0.0, 0.2, 0.4]) # To WU, Ex, BSM, Wo, Mov, Sl
    P[state_to_idx['Work'], :] = work_transitions / work_transitions.sum()

    # From 'Watch a Movie' (s4)
    P[state_to_idx['Watch a Movie'], state_to_idx['Work']] = 0.1
    P[state_to_idx['Watch a Movie'], state_to_idx['Watch a Movie']] = 0.7
    P[state_to_idx['Watch a Movie'], state_to_idx['Sleep']] = 0.2

    # From 'Sleep' (s5) - Terminal state
    P[state_to_idx['Sleep'], state_to_idx['Sleep']] = 1.0

    # 3. Solve the Bellman equation system: (I - gamma * P) * V = R
    I = np.identity(len(states))
    A = I - gamma * P
    V = np.linalg.solve(A, rewards)

    # 4. Output the results including the equations
    print("The state-value Bellman equation is V(s) = R(s) + γ * Σ [ P(s'|s) * V(s') ]")
    print("Here are the final equations for each state with the calculated values:\n")

    for i, state in enumerate(states):
        print(f"For state '{state}':")
        print(f"V({state}) = {V[i]:.8f}")
        
        # Build the right-hand side of the verification equation
        rhs_calc_str = f"{V[i]:.8f} = {rewards[i]:.1f} + {gamma} * ("
        
        sum_terms = []
        for j in range(len(states)):
            if P[i, j] > 0:
                sum_terms.append(f"{P[i, j]:.4f} * {V[j]:.8f}")
        
        rhs_calc_str += " + ".join(sum_terms) + ")"
        print(rhs_calc_str)
        print("-" * 30)

    # 5. Find and report the most promising state value
    most_promising_idx = np.argmax(V)
    most_promising_state = states[most_promising_idx]
    most_promising_value = V[most_promising_idx]

    print(f"\nThe most promising state is '{most_promising_state}'.")
    print(f"The state value of this most promising state is {most_promising_value:.8f}.")
    
    # 6. Return final answer in specified format
    print(f"\n<<<{most_promising_value:.8f}>>>")

if __name__ == "__main__":
    solve_working_life_problem()