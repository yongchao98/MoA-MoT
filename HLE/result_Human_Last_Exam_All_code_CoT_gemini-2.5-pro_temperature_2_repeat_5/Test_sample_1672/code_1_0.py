import numpy as np

def entropy(p):
    """Calculates the entropy of a probability distribution."""
    # Filter out probabilities equal to zero to avoid log(0) errors.
    p_positive = p[p > 1e-12]
    return -np.sum(p_positive * np.log(p_positive))

# 1. Define simulation parameters.
n_states = 20
n_iterations = 50
# `beta` is a learning rate that controls how aggressively the policy explores.
beta = 0.5

# 2. Initialize a non-uniform state distribution for the initial policy pi^0.
# This distribution is concentrated in the middle of the state space.
p_k = np.zeros(n_states)
p_k[n_states // 4 : 3 * n_states // 4] = 1.0
p_k = p_k / p_k.sum()

print(f"Demonstrating entropy evolution over {n_iterations} iterations for {n_states} states.")
print("-" * 60)

# 3. Track and print entropy at various iterations.
initial_entropy = entropy(p_k)
print(f"Entropy at k=0: {initial_entropy:.4f}")

# 4. Simulate the iterative process.
# At each step, we model the agent learning a new policy that flattens
# the state distribution. The update rule p_k ~ p_{k-1}^{(1-beta)} captures
# this effect of seeking novelty and making the distribution more uniform.
for k in range(1, n_iterations + 1):
    # Apply the update rule to get the new distribution.
    p_k = p_k**(1 - beta)
    p_k = p_k / p_k.sum()  # Re-normalize to ensure it's a valid probability distribution.

    # Print the entropy at key milestones to show the progression.
    if k in [1, 5, 10, 25, 50]:
        h_k = entropy(p_k)
        print(f"Entropy at k={k}: {h_k:.4f}")

# 5. Calculate the theoretical maximum entropy and compare with the final result.
max_entropy = np.log(n_states)
final_entropy = entropy(p_k)

print("-" * 60)
print(f"Final entropy after {n_iterations} iterations: {final_entropy:.4f}")
print(f"Theoretical maximum entropy (uniform distribution): {max_entropy:.4f}")
print("\nConclusion: As the number of iterations (k) increases, the policy pi^k")
print("induces a state distribution whose entropy converges to the maximum.")
