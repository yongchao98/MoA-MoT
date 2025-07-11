import numpy as np

def entropy(p):
    """Calculates the Shannon entropy of a probability distribution."""
    # Add a small epsilon to prevent log(0)
    p = p[p > 0]
    return -np.sum(p * np.log(p))

def softmax(x):
    """Compute softmax values for each sets of scores in x."""
    e_x = np.exp(x - np.max(x))
    return e_x / e_x.sum(axis=0)

def simulate_policy_updates():
    """
    Simulates the iterative policy update process to show that it maximizes state entropy.
    """
    # 1. Setup
    num_states = 5
    num_iterations = 200
    learning_rate = 0.5
    # Initialize with non-uniform logits, leading to a non-uniform policy
    logits = np.array([5.0, 1.0, 1.0, 1.0, 1.0])

    print("--- Simulation of Entropy Maximization ---")
    p_initial = softmax(logits)
    entropy_initial = entropy(p_initial)
    max_entropy = np.log(num_states)
    
    print(f"Number of states: {num_states}")
    print(f"Initial state distribution: {np.round(p_initial, 3)}")
    print(f"Initial entropy: {entropy_initial:.4f}")
    print(f"Maximum possible entropy: {max_entropy:.4f}\n")
    
    print("Running iterative updates...")

    # Store entropy history
    entropy_history = [entropy_initial]

    # 2. Iterative process
    p_k_minus_1 = p_initial
    current_logits = logits
    for k in range(1, num_iterations + 1):
        # Add a small epsilon for numerical stability
        epsilon = 1e-9
        
        # Reward is based on the previous policy's state distribution
        r_k = -np.log(p_k_minus_1 + epsilon)
        
        # Policy gradient update (simplified for this problem)
        # Gradient of J w.r.t logits is p_{k-1} * (r_k - E_{p_{k-1}}[r_k])
        # This pushes logits of low-prob states up and high-prob states down.
        expected_reward = np.sum(p_k_minus_1 * r_k)
        grad = p_k_minus_1 * (r_k - expected_reward)
        
        # Update logits
        current_logits += learning_rate * grad
        
        # The new policy/distribution
        p_k = softmax(current_logits)
        
        # Calculate and store new entropy
        entropy_k = entropy(p_k)
        entropy_history.append(entropy_k)
        
        # The new policy becomes the previous one for the next iteration
        p_k_minus_1 = p_k

        if k in [1, 5, 10, 50, 100, 200]:
             print(f"Iteration k={k:3d}: Current Entropy = {entropy_k:.4f}, Distribution = {np.round(p_k, 3)}")

    print("\n--- Conclusion ---")
    print("The simulation shows that as the number of iterations (k) increases, the state distribution becomes more uniform,")
    print("and its entropy approaches the maximum possible value.")
    print("This happens because the agent is continually incentivized to visit states that were previously infrequent.")
    print("\nThis means the policy that maximizes the state entropy is the limiting policy as k approaches infinity.")
    print("\nBased on this analysis, the correct option is A.")
    # The final equation is the convergence of the entropy to the maximum
    print(f"\nFinal State Equation:")
    final_p = np.round(p_k, 4)
    for i, p in enumerate(final_p):
        print(f"p(s_{i+1}) = {p}")
    print(f"Final Entropy H(s) = {entropy_history[-1]:.4f} ≈ log({num_states}) = {max_entropy:.4f}")


if __name__ == '__main__':
    simulate_policy_updates()
<<<A>>>