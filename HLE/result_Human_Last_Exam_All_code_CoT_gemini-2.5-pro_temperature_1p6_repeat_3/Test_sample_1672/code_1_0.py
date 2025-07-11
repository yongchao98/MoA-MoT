import numpy as np

def calculate_entropy_with_print(p, k):
    """Calculates the entropy of a probability distribution and prints the formula for the first step."""
    # Add a small epsilon to avoid log(0)
    p_safe = p + 1e-12
    entropy = -np.sum(p_safe * np.log(p_safe))

    if k == 0: # Print the equation for the first step as an example
        print("Example Entropy Calculation for k=0:")
        # Equation: H(p) = - Σ [p(s) * log(p(s))]
        equation_str = "H(p) = - ["
        terms = []
        for val in p:
            # Using f-string formatting to control the number of decimal places
            p_str = f"{val:.4f}"
            log_p_str = f"log({val:.4f})"
            terms.append(f"{p_str} * {log_p_str}")
        equation_str += " + ".join(terms) + "]"
        print(f"  {equation_str} = {entropy:.4f}")
        print("-" * 75)

    return entropy

# Main simulation logic
N_STATES = 5
max_entropy = np.log(N_STATES)

# Start with a skewed initial state distribution p^0(s) from an initial policy π^0.
p_k = np.array([0.5, 0.2, 0.15, 0.1, 0.05])
p_k = p_k / np.sum(p_k) # Normalize to ensure it's a valid probability distribution

print(f"This script simulates the evolution of the state distribution p_π^k(s).")
print(f"The reward at iteration k is r_k(s) = -log(p_π^(k-1)(s)).")
print("This encourages visiting less frequent states, which should increase state entropy over time.")
print(f"Maximum possible entropy for {N_STATES} states is log({N_STATES}) ≈ {max_entropy:.4f}.\n")

# A simplified model for the effect of the policy update.
# A new policy π^k will induce a distribution p^k that is 'flatter' than p^(k-1).
# We can model this with the update rule p_k ∝ (p_{k-1})^(1-β) for some β > 1,
# which captures the essence of shifting probability mass to low-probability states.
beta = 1.2
num_iterations = 15

for k in range(num_iterations):
    entropy = calculate_entropy_with_print(p_k, k)

    # Format the state distribution for printing
    dist_str = ", ".join([f"{val:.4f}" for val in p_k])
    print(f"k = {k:2d}: H(p) = {entropy:.4f}  |  p_π^k(s) = [{dist_str}]")

    # The rewards r_{k+1} = -log(p_k) are used to train the next policy π^{k+1}.
    # The resulting distribution p_{k+1} is then calculated for the next step.
    p_unnormalized = p_k**(1 - beta)
    p_k = p_unnormalized / np.sum(p_unnormalized)

print("\n" + "="*75)
print("Conclusion:")
print("The simulation shows that as k increases, the state distribution p_π^k(s) converges")
print(f"towards a uniform distribution, and its entropy H(p) approaches the maximum value of {max_entropy:.4f}.")
print("This confirms that the policy which maximizes state entropy is the limiting policy obtained")
print("after a large number of iterations.")