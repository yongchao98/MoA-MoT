import numpy as np

def solve_mrp():
    """
    Solves the state-value for the given Markov Reward Process.
    """
    # 1. Define the components of the Markov Reward Process
    states = ["Wake Up", "Exercise", "Browse Social Media", "Work", "Watch a Movie", "Sleep"]
    n_states = len(states)

    # Reward vector R
    R = np.array([1, 3, -2, 2, 1, 0], dtype=float)

    # Discount factor gamma
    gamma = 0.2

    # 2. Define the transition probability matrix P
    # NOTE: Probabilities from "Work" in the graph sum to 1.1 (0.1+0.1+0.3+0.2+0.4).
    # We normalize them to sum to 1.0.
    sum_prob_work = 0.1 + 0.1 + 0.3 + 0.2 + 0.4

    P = np.array([
        # To: WakeUp, Exercise, B.Social, Work, W.Movie, Sleep
        [0.0,   0.2,    0.5,    0.1,    0.1,    0.1],      # From Wake Up
        [0.0,   0.0,    0.2,    0.5,    0.0,    0.3],      # From Exercise
        [0.0,   0.4,    0.0,    0.6,    0.0,    0.0],      # From Browse Social Media
        [0.1/sum_prob_work, 0.1/sum_prob_work, 0.3/sum_prob_work, 0.0, 0.2/sum_prob_work, 0.4/sum_prob_work], # From Work (Normalized)
        [0.0,   0.0,    0.0,    0.1,    0.7,    0.2],      # From Watch a Movie
        [0.0,   0.0,    0.0,    0.0,    0.0,    0.0]       # From Sleep (Terminal state)
    ], dtype=float)

    # 3. Solve the Bellman equation: V = (I - gamma * P)^-1 @ R
    I = np.identity(n_states)
    try:
        V = np.linalg.solve(I - gamma * P, R)
    except np.linalg.LinAlgError:
        print("Error: The matrix is singular, and the system cannot be solved directly.")
        return

    # 4. Print the detailed equations and the solution
    print("The state-value Bellman equation is V(s) = R(s) + γ * Σ [P(s'|s) * V(s')].")
    print("This creates a system of linear equations which we solve for V.\n")

    print("--- Equations ---")
    for i in range(n_states):
        if states[i] == 'Sleep':
            print(f"V(Sleep) = 0")
            continue
        equation = f"V({states[i]}) = {R[i]:.1f} + {gamma} * ("
        terms = [f"{P[i, j]:.4f}*V({states[j]})" for j in range(n_states) if P[i, j] > 0]
        equation += " + ".join(terms) + ")"
        print(equation)

    print("\n--- Solved State Values ---")
    for i in range(n_states):
        print(f"V({states[i]:<20}) = {V[i]:.8f}")
    
    print("\n--- Verification of a Sample Equation (for 'Work') ---")
    i = states.index('Work')
    sum_term = np.dot(P[i, :], V)
    verification_eq = f"{V[i]:.8f} = {R[i]:.1f} + {gamma} * ("
    terms = []
    for j in range(n_states):
        if P[i,j] > 0:
            terms.append(f"{P[i, j]:.4f} * {V[j]:.8f}")
    verification_eq += " + ".join(terms)
    verification_eq += f")\n{V[i]:.8f} = {R[i]:.1f} + {gamma} * {sum_term:.8f}"
    print(verification_eq)
    print(f"{V[i]:.8f} = {R[i] + gamma * sum_term:.8f}")

    # Find and report the most promising state
    most_promising_idx = np.argmax(V)
    most_promising_state = states[most_promising_idx]
    most_promising_value = V[most_promising_idx]

    print(f"\nThe most promising state is '{most_promising_state}'.")
    print(f"The value of this state is {most_promising_value:.8f}.")
    
    # Return final answer in specified format
    print(f"\n<<<{most_promising_value:.8f}>>>")

if __name__ == "__main__":
    solve_mrp()