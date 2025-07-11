import numpy as np

# This script calculates the state values for a given Markov Reward Process
# using the analytical solution to the Bellman equation.

# --- 1. Problem Definition ---

# State names mapping for clarity in the output
states = {
    0: "Wake Up",
    1: "Exercise",
    2: "Browse Social Media",
    3: "Work",
    4: "Watch a Movie",
    5: "Sleep"
}
num_states = len(states)

# Discount factor
gamma = 0.2

# Rewards for each state (R)
# R = [R(Wake Up), R(Exercise), R(BSM), R(Work), R(Movie), R(Sleep)]
R = np.array([1, 3, -2, 2, 1, 0], dtype=float)

# Transition Probability Matrix (P), where P[i, j] is the prob of moving from state i to j.
# Note: The sum of probabilities for outgoing transitions from 'Work' in the graph is 1.1.
# We normalize these probabilities so they sum to 1.
work_probs_raw = np.array([0.1, 0.1, 0.3, 0.0, 0.2, 0.4]) # To WU, Ex, BSM, Work, Movie, Sleep
work_probs_normalized = work_probs_raw / work_probs_raw.sum()

P = np.array([
    # WU      Ex      BSM     Work     Movie   Sleep
    [0.0,     0.2,    0.5,    0.1,     0.1,    0.1],    # from Wake Up
    [0.0,     0.0,    0.2,    0.5,     0.0,    0.3],    # from Exercise
    [0.0,     0.4,    0.0,    0.6,     0.0,    0.0],    # from Browse Social Media
    [work_probs_normalized[0], work_probs_normalized[1], work_probs_normalized[2], 0.0, work_probs_normalized[4], work_probs_normalized[5]], # from Work (normalized)
    [0.0,     0.0,    0.0,    0.1,     0.7,    0.2],    # from Watch a Movie
    [0.0,     0.0,    0.0,    0.0,     0.0,    1.0]     # from Sleep (terminal state)
], dtype=float)


# --- 2. Solving the Bellman Equation ---

# The Bellman equation in matrix form is: V = R + gamma * P * V
# We solve this system of linear equations: (I - gamma * P) * V = R
# V = (I - gamma * P)^-1 * R
I = np.identity(num_states)
V = np.linalg.inv(I - gamma * P) @ R


# --- 3. Displaying the Results ---

print("The solved Bellman equations for each state are:\n")

for i in range(num_states):
    state_name = states[i]
    reward = R[i]
    
    # Build the string for the sum part of the equation: Σ P(s'|s) * V(s')
    sum_part_str = " + ".join([f"{P[i, j]:.9f} * {V[j]:.8f}" for j in range(num_states) if P[i, j] > 0])
    
    # Print the equation with values plugged in
    print(f"V({state_name}) = {reward:.1f} + {gamma} * ({sum_part_str}) = {V[i]:.8f}\n")


# --- 4. Final Answer ---

# Find the state with the highest value
max_value = np.max(V)
max_index = np.argmax(V)
most_promising_state_name = states[max_index]

print("-" * 60)
print(f"The most promising state is '{most_promising_state_name}'.")
print(f"The state value of the most promising state is: {max_value:.8f}")
print("-" * 60)