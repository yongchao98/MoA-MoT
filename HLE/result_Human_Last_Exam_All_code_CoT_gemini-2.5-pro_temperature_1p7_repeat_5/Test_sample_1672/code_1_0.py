import numpy as np

def entropy(p):
    """Calculates the entropy of a probability distribution."""
    # Add a small epsilon to prevent log(0)
    epsilon = 1e-12
    return -np.sum(p * np.log(p + epsilon))

def run_simulation():
    """
    Simulates the iterative process of entropy maximization and prints the results.
    """
    # --- Simulation Parameters ---
    N_STATES = 20      # Number of states in the environment
    ITERATIONS = 200   # Number of iterations to run
    ALPHA = 0.5        # Learning rate for updating the distribution
    EPSILON = 1e-9     # Small value to avoid division by zero

    print("--- Entropy Maximization Simulation ---")
    print(f"Number of states: {N_STATES}")
    
    # 1. Initialize with a highly skewed (low entropy) distribution
    p = np.zeros(N_STATES)
    p[0] = 1.0
    
    initial_entropy = entropy(p)
    print(f"\nInitial Entropy (at k=0): {initial_entropy:.4f}")

    # 2. The maximum possible entropy is for a uniform distribution
    max_entropy = np.log(N_STATES)
    print(f"Maximum possible Entropy: {max_entropy:.4f}\n")
    
    print("--- Running Iterations ---")
    # 3. Iteratively update the distribution
    iterations_to_print = [1, 2, 5, 10, 20, 50, 100, ITERATIONS]
    for k in range(1, ITERATIONS + 1):
        # The reward is r = -log(p_old). A policy maximizing this will visit states
        # with low p_old. We can model the resulting *target* distribution
        # as being proportional to 1/p_old.
        target_dist_unnormalized = 1 / (p + EPSILON)
        target_dist = target_dist_unnormalized / np.sum(target_dist_unnormalized)

        # Update the current distribution to move towards the target (like a policy update)
        p = (1 - ALPHA) * p + ALPHA * target_dist
        
        # Print entropy at specific iterations
        if k in iterations_to_print:
            current_entropy = entropy(p)
            print(f"Entropy at iteration k={k:<4}: {current_entropy:.4f}")

    print("\n--- Conclusion ---")
    print("The simulation shows that with each iteration (k), the policy is updated")
    print("to explore less-visited states, making the state distribution more uniform.")
    print("As a result, the entropy of the state distribution consistently increases")
    print("and approaches the theoretical maximum entropy.")
    print("\nTherefore, the policy that maximizes the entropy is the one obtained at the limit")
    print("of this process, as k approaches infinity.")


if __name__ == '__main__':
    run_simulation()
    # The final answer is A, corresponding to the limiting policy.
    # Appending the answer in the required format.
    print("\n<<<A>>>")
