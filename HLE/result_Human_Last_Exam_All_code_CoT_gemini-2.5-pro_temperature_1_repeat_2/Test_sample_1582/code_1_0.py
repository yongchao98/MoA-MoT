import random

# This script simulates a Markov chain that satisfies the conditions of the problem
# to demonstrate that it is not positive recurrent.
# The chosen example is a biased random walk on the non-negative integers Z+.
# - State space Sigma = {0, 1, 2, ...}
# - Transition probabilities: p(i, i+1) = p, p(i, i-1) = 1-p for i > 0, and p(0, 1) = 1.
# - We choose a bias p > 0.5, which creates a drift away from the origin.
# This chain is irreducible.

# Parameters for the simulation
# Probability of moving to the right (for i > 0)
p_right = 0.75
# Probability of moving to the left (for i > 0)
p_left = 1 - p_right

# The problem conditions are met with:
# - Finite set A = {0}.
# - Function f(x) = (p_right / p_left)**x. For p_right > 0.5, this function tends to infinity.
# - For x > 0 (i.e., x not in A), the drift condition E[f(X_n+1)|X_n=x] - f(x) >= 0 holds.

# We start the walk from a state outside A.
start_state = 1

# Simulation settings
num_walks = 100000
max_steps_per_walk = 1000

# --- Simulation ---
returned_to_zero = 0

for _ in range(num_walks):
    current_state = start_state
    for _ in range(max_steps_per_walk):
        # If the walker returns to state 0 (the set A), we count it and stop this walk.
        if current_state == 0:
            returned_to_zero += 1
            break
        
        # Transition logic for the biased random walk
        if random.random() < p_right:
            current_state += 1
        else:
            current_state -= 1
    # If the loop finishes, the walker did not return to 0 within the step limit.
    # For a transient chain with drift to +inf, this is the likely outcome.

# --- Output Results ---
estimated_prob = returned_to_zero / num_walks
# For a biased random walk starting at 1, the theoretical probability of ever hitting 0
# is min(1, (1-p)/p).
theoretical_prob = min(1.0, p_left / p_right)

print("--- Markov Chain Simulation: Biased Random Walk ---")
print("This simulation models a chain satisfying the problem's conditions.")
print("A biased random walk with drift to infinity is transient, thus not positive recurrent.")
print(f"We simulate walks from state {start_state} and check if they return to A={{0}}.")
print("\n--- Simulation Results ---")
print(f"Number of simulated walks: {num_walks}")
print(f"Number of walks that returned to state 0: {returned_to_zero}")
print(f"Estimated probability of returning to 0: {estimated_prob:.4f}")

print("\n--- Theoretical Calculation ---")
print(f"The theoretical probability of a walk starting at {start_state} ever returning to 0 is (1-p)/p.")
print(f"Equation: P(return) = (1 - {p_right}) / {p_right}")
print(f"Final equation result: {p_left} / {p_right} = {theoretical_prob:.4f}")

print("\n--- Final Conclusion ---")
print("The simulation shows the probability of return is < 1, supporting the conclusion")
print("that the chain is transient and therefore not positive recurrent.")
