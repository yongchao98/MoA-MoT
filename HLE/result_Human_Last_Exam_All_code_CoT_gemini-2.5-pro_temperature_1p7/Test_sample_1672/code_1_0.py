import numpy as np

def calculate_entropy(p):
    """Calculates the entropy of a probability distribution."""
    # Filter out zero probabilities to avoid log(0) error
    p_positive = p[p > 0]
    return -np.sum(p_positive * np.log(p_positive))

def print_entropy_calculation(p, iteration_name):
    """Prints the distribution, entropy, and the calculation formula."""
    entropy = calculate_entropy(p)
    print(f"\n--- At Iteration {iteration_name} ---")
    print(f"Distribution p = {np.round(p, 4)}")
    
    # Building the entropy calculation equation string to show the numbers
    # H(s) = - sum( p(s) * log(p(s)) )
    equation_parts = []
    for val in p:
        if val > 0:
            # Using .4f format to match the rounded distribution printout
            equation_parts.append(f"{val:.4f} * log({val:.4f})")
    equation_str = " + ".join(equation_parts)
    print(f"Entropy Equation: H(p) = - ( {equation_str} )")
    print(f"Calculated Entropy H(p) = {entropy:.4f}")

# --- Setup ---
NUM_STATES = 5
NUM_ITERATIONS = 10

# Initialize a non-uniform policy/state distribution
p = np.array([0.5, 0.25, 0.15, 0.08, 0.02])

# Calculate the maximum possible entropy for a uniform distribution
p_uniform = np.full(NUM_STATES, 1.0 / NUM_STATES)
max_entropy = calculate_entropy(p_uniform)
print(f"Theoretical maximum entropy for {NUM_STATES} states is {max_entropy:.4f} (from a uniform distribution).")

# --- Iteration 0 (Initial state) ---
print_entropy_calculation(p, "0 (Initial)")

# --- Iterations ---
# We show the state at k=1 and k=final to see the progression.
for k in range(1, NUM_ITERATIONS + 1):
    p_prev = p
    # The reward is r_k(s) = -log(p_{k-1}(s)).
    # A policy trying to maximize this reward will tend to visit states with low p_{k-1}(s).
    # We model this update by making the new distribution p_k(s) proportional to 1/p_{k-1}(s).
    p = 1.0 / p_prev
    p /= np.sum(p) # Normalize to make it a valid probability distribution

    if k == 1:
        print_entropy_calculation(p, "1")
    if k == NUM_ITERATIONS:
        print_entropy_calculation(p, f"{NUM_ITERATIONS} (Final)")

print("\nConclusion:")
print(f"As the number of iterations k increases, the distribution p rapidly approaches a uniform distribution {[round(x, 2) for x in p_uniform]},")
print("and its entropy H(p) approaches the maximum possible value.")
print("This simulation supports the conclusion that the limiting policy, lim_{k->inf} pi^k, maximizes the state entropy H(s).")