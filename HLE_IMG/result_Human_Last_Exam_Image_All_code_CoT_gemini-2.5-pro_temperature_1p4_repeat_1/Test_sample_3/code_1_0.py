import numpy as np

def solve_state_values():
    """
    Calculates the state values for the given Markov Reward Process using the Bellman equation.
    """
    # Define states, rewards, and the discount factor
    states = ["Wake Up", "Exercise", "Browse Social Media", "Work", "Watch a Movie", "Sleep"]
    state_vars = ["V(Wake Up)", "V(Exercise)", "V(Browse)", "V(Work)", "V(Movie)", "V(Sleep)"]
    rewards = np.array([1, 3, -2, 2, 1, 0])
    gamma = 0.2

    # Define the transition probability matrix P from the graph.
    # Note: The outgoing probabilities from 'Work' in the graph sum to 1.1.
    # (0.1+0.3+0.1+0.2+0.4=1.1). We assume the transition Work->Sleep is 0.3
    # instead of 0.4 to make the probabilities sum to 1.0.
    P = np.array([
        #  WU,  EX, BSM,   W,  WM,   S
        [ 0.0, 0.2, 0.5, 0.1, 0.1, 0.1],  # Wake Up
        [ 0.0, 0.0, 0.2, 0.5, 0.0, 0.3],  # Exercise
        [ 0.0, 0.4, 0.0, 0.6, 0.0, 0.0],  # Browse Social Media
        [ 0.1, 0.1, 0.3, 0.0, 0.2, 0.3],  # Work (corrected)
        [ 0.0, 0.0, 0.0, 0.1, 0.7, 0.2],  # Watch a Movie
        [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]   # Sleep (Terminal state)
    ])

    print("Bellman equations V(s) = R(s) + gamma * sum(P(s'|s) * V(s')):\n")
    for i in range(len(states)):
        equation = f"{state_vars[i]} = {rewards[i]} + {gamma} * ("
        terms = []
        for j in range(len(states)):
            if P[i, j] > 0:
                terms.append(f"{P[i, j]} * {state_vars[j]}")
        if not terms:
            equation += "0"
        else:
            equation += " + ".join(terms)
        equation += ")"
        print(equation)

    # Solve the Bellman equation system: (I - gamma * P) * V = R
    I = np.identity(len(states))
    A = I - gamma * P
    
    # Use np.linalg.solve which is more numerically stable than finding the inverse
    V = np.linalg.solve(A, rewards)

    print("\nCalculated State Values:")
    for i in range(len(states)):
        print(f"{state_vars[i]:<20} = {V[i]:.8f}")

    # Find the most promising state
    max_value = np.max(V)
    most_promising_state_index = np.argmax(V)
    most_promising_state_name = states[most_promising_state_index]

    print(f"\nThe most promising state is '{most_promising_state_name}' with a value of {max_value:.8f}.")

    print(f"\nFinal Answer: The state value of the most promising state is {max_value:.8f}")

if __name__ == '__main__':
    solve_state_values()
    # To directly output the answer for the platform
    # This part is for programmatic extraction of the answer
    rewards = np.array([1, 3, -2, 2, 1, 0])
    gamma = 0.2
    P = np.array([
        [0.0, 0.2, 0.5, 0.1, 0.1, 0.1],
        [0.0, 0.0, 0.2, 0.5, 0.0, 0.3],
        [0.0, 0.4, 0.0, 0.6, 0.0, 0.0],
        [0.1, 0.1, 0.3, 0.0, 0.2, 0.3],
        [0.0, 0.0, 0.0, 0.1, 0.7, 0.2],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    ])
    I = np.identity(6)
    A = I - gamma * P
    V = np.linalg.solve(A, rewards)
    max_value = np.max(V)
    print(f"<<<{max_value:.8f}>>>")