import math
import numpy as np

def calculate_and_print_entropy(p, name="p"):
    """Calculates and prints the entropy equation and result for a given distribution."""
    if not np.isclose(sum(p), 1.0):
        raise ValueError("Probabilities must sum to 1.")
    
    entropy = -sum([x * math.log(x) for x in p if x > 0])
    
    # Building the equation string
    equation_parts = []
    for x in p:
        if x > 0:
            equation_parts.append(f"{x:.3f}*ln({x:.3f})")
    
    equation_str = f"H({name}) = - ( " + " + ".join(equation_parts) + " )"
    
    print(f"Let's assume the state distribution for a policy {name} is {np.round(p, 3)}.")
    print(f"The entropy calculation is:")
    print(f"{equation_str} = {entropy:.4f}")
    print("-" * 20)

# --- Main Demonstration ---

# Step 0: Initial, non-uniform distribution from policy π^0
p0 = np.array([0.1, 0.8, 0.1])
calculate_and_print_entropy(p0, "p_π^0")

# Step 1: Policy π^1 is rewarded for visiting states s1 and s3, making the distribution more uniform.
p1 = np.array([0.3, 0.4, 0.3])
calculate_and_print_entropy(p1, "p_π^1")

# The limit: A uniform distribution that would be induced by lim (k->inf) π^k
p_limit = np.array([1/3, 1/3, 1/3])
calculate_and_print_entropy(p_limit, "p_limit")

print("As the example shows, H(p_π^1) > H(p_π^0), and the entropy is maximized at the uniform distribution (p_limit).")
print("This demonstrates that the iterative process drives the policy towards maximizing state entropy, which is achieved in the limit.")
print("\nTherefore, the policy that maximizes the entropy H(s) is lim_{k->inf} π^k.")