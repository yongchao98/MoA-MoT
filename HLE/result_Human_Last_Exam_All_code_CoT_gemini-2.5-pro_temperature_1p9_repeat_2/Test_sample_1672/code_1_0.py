import numpy as np

def entropy(p):
    """Calculates the entropy of a probability distribution."""
    # Filter out zero probabilities to avoid log(0)
    p = p[p > 1e-12]
    return -np.sum(p * np.log(p))

# --- Simulation Parameters ---
NUM_STATES = 10
NUM_ITERATIONS = 20
# Learning rate 'alpha' for the update. Controls how fast p converges.
ALPHA = 0.5

# --- Initialization ---
# Start with a skewed, non-uniform distribution for p at k=0.
# This represents the state distribution p_{\pi^0}(s).
p = np.array([0.5, 0.2, 0.1, 0.05, 0.05, 0.02, 0.02, 0.02, 0.02, 0.02])
p /= np.sum(p)

print(f"Demonstrating convergence to the maximum entropy policy.\n")
print(f"Starting with a skewed distribution p_0 with {NUM_STATES} states:")
print(f"p_0 = {np.round(p, 3)}")
print(f"Initial Entropy H(p_0) = {entropy(p):.4f}\n" + "-"*30)

# --- Iterative Process ---
for k in range(1, NUM_ITERATIONS + 1):
    # The current distribution `p` is p_{k-1}
    
    # 1. Calculate the reward r_k(s) = -log(p_{k-1}(s))
    # We clip p to avoid numerical issues with very small values.
    p_clipped = np.clip(p, 1e-9, 1.0)
    r = -np.log(p_clipped)
    
    # 2. The new policy pi^k aims to maximize this reward. Its resulting
    #    distribution p_k will shift mass towards states with high reward.
    #    We model this by defining a "target" distribution proportional to exp(r),
    #    which is equivalent to 1/p.
    target_dist_unnormalized = 1 / p_clipped
    target_dist = target_dist_unnormalized / np.sum(target_dist_unnormalized)

    # 3. The new distribution p_k is a mix of the old p_{k-1} and the target.
    #    This simulates a gradual policy update.
    p = (1 - ALPHA) * p + ALPHA * target_dist

    # --- Reporting for the current iteration ---
    print(f"Iteration k={k}")
    print(f"Entropy H(p_{k}) = {entropy(p):.4f}")
    if k % 5 == 0 or k == NUM_ITERATIONS:
      print(f"p_{k} = {np.round(p, 3)}")
      print("-"*30)


# --- Final Result ---
final_p = p
final_entropy = entropy(final_p)

# The maximum possible entropy for this many states is log(NUM_STATES)
max_entropy = np.log(NUM_STATES)

print(f"\nAfter {NUM_ITERATIONS} iterations, the entropy is {final_entropy:.4f}.")
print(f"The theoretical maximum entropy for {NUM_STATES} states is {max_entropy:.4f}.")
print("The simulation shows the entropy converging towards this maximum.")
print("This corresponds to the policy lim_{k->inf} pi^k, which induces a uniform state distribution.")

print("\n--- Final Entropy Calculation ---")
print("The final distribution is p =", np.round(final_p, 3))
print("The entropy is calculated as: H(p) = - sum( p[s] * log(p[s]) ) for all states s.")
# Build and print the full calculation string
equation_str = "H(p_final) = -(" + " + ".join([f"{val:.3f}*log({val:.3f})" for val in final_p]) + ")"
print("For our final distribution, the calculation is:")
print(equation_str)
print(f"Result: H(p_final) = {final_entropy:.4f}")
