import numpy as np

def calculate_entropy(p, p_name):
    """Calculates and prints the entropy of a probability distribution."""
    # Filter out zero probabilities to avoid log(0) and ensure the sum is 1.
    p = p[p > 0]
    p = p / np.sum(p)
    
    entropy_val = -np.sum(p * np.log(p))
    
    # Building the equation string
    equation_parts = [f"({pi:.2f} * log({pi:.2f}))" for pi in p]
    equation_str = " + ".join(equation_parts)
    
    print(f"For distribution {p_name}:")
    print(f"  p = {np.round(p, 4)}")
    print(f"  Entropy H = -[ {equation_str} ] = {entropy_val:.4f}\n")
    return entropy_val

# Let's consider an environment with 5 states.

# Iteration 0: The initial policy pi^0 induces a skewed state distribution.
p0 = np.array([0.8, 0.1, 0.05, 0.03, 0.02])
calculate_entropy(p0, "p_pi^0")

# Iteration 1: The policy pi^1 is rewarded for visiting novel states,
# leading to a more uniform distribution.
p1 = np.array([0.4, 0.25, 0.2, 0.1, 0.05])
calculate_entropy(p1, "p_pi^1")

# Iteration k: After several iterations, the distribution is even more uniform.
pk = np.array([0.22, 0.21, 0.20, 0.19, 0.18])
calculate_entropy(pk, "p_pi^k")

# Limit k -> infinity: The process converges to a policy that induces a
# perfectly uniform distribution, which maximizes entropy.
p_inf = np.array([0.2, 0.2, 0.2, 0.2, 0.2])
calculate_entropy(p_inf, "p_pi^inf")

print("As shown, as the state distribution becomes more uniform, its entropy increases.")
print("The iterative process is designed to make the distribution more uniform at each step.")
print("Therefore, the maximum entropy is achieved in the limit as k -> infinity.")
