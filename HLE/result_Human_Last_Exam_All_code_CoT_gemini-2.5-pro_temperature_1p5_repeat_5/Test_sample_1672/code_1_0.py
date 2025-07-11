import numpy as np

def calculate_entropy(p):
    """Calculates the entropy (in bits) of a probability distribution."""
    # Filter out zero probabilities to avoid log(0) which is undefined.
    # In practice, p(s) > 0 for any reachable state.
    p = p[p > 0]
    return -np.sum(p * np.log2(p))

# --- Setup ---
# Imagine an environment with 5 states. We'll compare the entropy
# of different state visitation distributions.

# Distribution 1: Skewed distribution (like from an initial policy pi^0)
# The agent starts by mostly visiting the first state.
p_initial = np.array([0.8, 0.1, 0.05, 0.03, 0.02])

# Distribution 2: An intermediate distribution (like from a policy pi^k)
# The agent has learned to explore other states to get the intrinsic reward.
p_intermediate = np.array([0.4, 0.25, 0.15, 0.1, 0.1])

# Distribution 3: A uniform distribution (the limit case)
# The converged policy lim_{k->inf} pi^k induces a uniform state distribution.
p_uniform = np.array([0.2, 0.2, 0.2, 0.2, 0.2])

# --- Calculations & Output ---
entropy_initial = calculate_entropy(p_initial)
entropy_intermediate = calculate_entropy(p_intermediate)
entropy_uniform = calculate_entropy(p_uniform)

print("This script demonstrates that state entropy H(s) is maximized by a uniform state distribution.")
print("The learning process described uses rewards that encourage visiting rare states,")
print("driving the state distribution towards uniformity as the number of iterations k increases.\n")

print("-" * 50)
print(f"1. A skewed distribution (representing an early policy like pi^0):")
print(f"   p(s) = {p_initial}")
print(f"   Entropy H(s) = - (", " + ".join([f"{p:.2f}*log2({p:.2f})" for p in p_initial]), f") = {entropy_initial:.4f} bits")
print("-" * 50)

print(f"2. An intermediate distribution (representing a later policy pi^k):")
print(f"   p(s) = {p_intermediate}")
print(f"   Entropy H(s) = - (", " + ".join([f"{p:.2f}*log2({p:.2f})" for p in p_intermediate]), f") = {entropy_intermediate:.4f} bits")
print("-" * 50)

print(f"3. A uniform distribution (representing the converged policy lim_{k->inf} pi^k):")
print(f"   p(s) = {p_uniform}")
print(f"   Entropy H(s) = - (", " + ".join([f"{p:.2f}*log2({p:.2f})" for p in p_uniform]), f") = {entropy_uniform:.4f} bits")
print("-" * 50)

print("\nAs shown, the entropy is highest for the uniform distribution.")
print("Therefore, the policy that maximizes the state entropy H(s) is the one the process converges to.")
