import numpy as np

# Define the states for labeling the output
states = ["Wake Up", "Exercise", "Browse Social Media", "Work", "Watch a Movie", "Sleep"]

# Define the reward vector for each state
R = np.array([1, 3, -2, 2, 1, 0])

# Define the discount factor
gamma = 0.2

# Define the transition probability matrix P(s'|s)
# Rows are current states, columns are next states.
P = np.array([
    #       WU,  Ex,  BSM, Work,Movie,Sleep
    [0.0, 0.2, 0.5, 0.1, 0.1, 0.1],  # From Wake Up
    [0.0, 0.0, 0.2, 0.5, 0.0, 0.3],  # From Exercise
    [0.0, 0.4, 0.0, 0.6, 0.0, 0.0],  # From Browse Social Media
    [0.0, 0.1, 0.3, 0.0, 0.2, 0.4],  # From Work (assuming 0.4 to Sleep)
    [0.0, 0.0, 0.0, 0.1, 0.7, 0.2],  # From Watch a Movie
    [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]   # From Sleep (Terminal)
])

# Solve the Bellman equation: V = (I - gamma*P)^-1 * (P*R)
I = np.identity(len(states))
A = I - gamma * P
B = P @ R

# Solve the system of linear equations A*V = B
V = np.linalg.solve(A, B)

# Find the most promising state (the one with the highest value)
max_value = np.max(V)
max_idx = np.argmax(V)
most_promising_state = states[max_idx]

# Print all state values for context
print("Calculated State Values (V):")
for i, state in enumerate(states):
    print(f"  V({state:20s}) = {V[i]:.8f}")
print("-" * 40)

# Print the Bellman equation for the most promising state
print(f"The most promising state is '{most_promising_state}' with a value of {max_value:.8f}.\n")
print(f"The Bellman equation for V({most_promising_state}) is:")
print(f"V(s) = Σ P(s'|s) * [R(s') + γ * V(s')]")

# Build and print the equation with all numbers plugged in
equation_str = f"{V[max_idx]:.8f} = "
equation_terms = []
for next_state_idx, prob in enumerate(P[max_idx]):
    if prob > 0:
        term = f"{prob:.1f} * ({R[next_state_idx]} + {gamma} * {V[next_state_idx]:.8f})"
        equation_terms.append(term)

equation_str += " + ".join(equation_terms)
print(equation_str)

# Finally, output the most promising state's value as the answer
print(f"\nThe state value of the most promising state is {max_value:.8f}.")
