import numpy as np

def calculate_entropy(p):
    """Calculates the entropy of a probability distribution."""
    # Filter out zero probabilities to avoid log(0)
    p_nz = p[p > 0]
    return -np.sum(p_nz * np.log(p_nz))

def main():
    """
    Simulates the evolution of a state distribution driven by entropy-seeking rewards.
    """
    # Define environment and simulation parameters
    num_states = 10
    num_iterations = 40
    # The parameter 'alpha' controls the rate of change. This simplified update rule
    # p_k ~ (p_{k-1})^(1-alpha) captures the essence of moving towards states with
    # high reward r_k = -log(p_{k-1}).
    alpha = 0.5
    # Small constant to prevent log(0) errors
    epsilon = 1e-12

    print("This simulation shows how a state distribution evolves over iterations.")
    print("The update rule encourages visiting less probable states, which drives")
    print("the distribution towards uniformity, thereby maximizing its entropy.")
    print("-" * 70)

    # Start with a highly skewed distribution for k=0 (e.g., policy pi^0 always stays at state 0)
    p = np.zeros(num_states)
    p[0] = 1.0

    entropy_at_0 = calculate_entropy(p)
    print(f"Iteration k=0: Entropy = {entropy_at_0:.4f}")

    # Run the iterative process
    for k in range(1, num_iterations + 1):
        # The reward for iteration k is r_k(s) = -log(p_{k-1}(s)).
        # The new policy pi^k will result in a new distribution p_k(s).
        # We model this change by making the new probability proportional to a power
        # of the old probability, where the exponent is less than 1.
        # This gives higher relative weight to states with lower initial probability.
        
        # Add epsilon to handle states where p=0.
        p_prev_safe = p + epsilon

        # Update the distribution
        p_new = p_prev_safe**(1 - alpha)
        
        # Normalize to ensure it remains a valid probability distribution
        p = p_new / np.sum(p_new)
        
        # Calculate and print entropy for the new distribution
        entropy_k = calculate_entropy(p)
        if k % 5 == 0 or k == 1:
            print(f"Iteration k={k:<2}: Entropy = {entropy_k:.4f}")

    # Calculate the theoretical maximum entropy for a uniform distribution
    max_entropy = np.log(num_states)
    
    print("-" * 70)
    print(f"Final Entropy at k={num_iterations}: {entropy_k:.4f}")
    print(f"Theoretical Maximum Entropy: {max_entropy:.4f}")
    print("\nAs k increases, the entropy converges to its maximum possible value.")
    print("This confirms that the policy that maximizes the state entropy H(s) is the limit of pi^k as k tends to infinity.")

if __name__ == '__main__':
    main()
