import numpy as np

def calculate_entropy(p):
    """Calculates the entropy of a probability distribution."""
    # Add a small epsilon to prevent log(0)
    p = p + 1e-9
    return -np.sum(p * np.log(p))

def run_simulation():
    """
    Simulates the iterative process of entropy maximization.
    """
    # 1. Define simulation parameters
    num_states = 5
    num_iterations = 10
    # The 'alpha' parameter controls the learning rate. For convergence, 0 < alpha < 2.
    # A value of 1 would cause oscillation in this model. We choose 0.8 for stable convergence.
    alpha = 0.8

    # 2. Create a non-uniform initial state distribution for the policy pi^0
    # This represents the state distribution p_{pi^0}(s)
    p = np.array([0.7, 0.15, 0.05, 0.05, 0.05])
    
    print(f"Starting simulation with {num_states} states.\n")
    print(f"Maximum possible entropy (for a uniform distribution) is: {np.log(num_states):.4f}\n")
    
    # 3. Run the iterative process
    for k in range(num_iterations):
        # Calculate the entropy of the current state distribution p_k
        entropy = calculate_entropy(p)
        
        # Policy pi^k induces distribution p_k. Let's analyze it.
        # This corresponds to the distributions from policies pi^0, pi^1, ..., pi^(k)
        print(f"--- Iteration k = {k} (Policy pi^{k}) ---")
        print(f"State Distribution p_{k}(s): {np.round(p, 4)}")
        print(f"Entropy H(s): {entropy:.4f}")
        print("-" * (35 + len(str(k))))
        
        # Update the policy for the next iteration k+1.
        # The reward for the new policy pi^{k+1} is r_{k+1}(s) = -log(p_k(s)).
        # The new policy pi^{k+1} will try to visit states to maximize this reward.
        # We model this by making the new state distribution p_{k+1}(s)
        # proportional to exp(reward), which simplifies to (p_k(s))^{-alpha}.
        
        # Add a small constant to p to avoid division by zero if a probability is 0
        p_safe = p + 1e-9
        
        # Unnormalized new distribution
        p_unnormalized = p_safe**(-alpha)
        
        # Normalize to get the new probability distribution p_{k+1}
        p = p_unnormalized / np.sum(p_unnormalized)

    print("\nConclusion:")
    print("As the number of iterations 'k' increases, the entropy of the state distribution")
    print("induced by the policy pi^k steadily increases and approaches the maximum possible value.")
    print("This demonstrates that the policy that maximizes the entropy is the one found in the limit")
    print("as k -> infinity.")

if __name__ == '__main__':
    run_simulation()
<<<A>>>