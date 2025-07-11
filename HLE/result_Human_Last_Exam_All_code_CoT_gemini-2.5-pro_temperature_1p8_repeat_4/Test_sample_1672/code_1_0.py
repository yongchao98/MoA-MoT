import numpy as np

def calculate_entropy(p):
    """Calculates the entropy of a probability distribution in bits."""
    # Filter out zero probabilities to avoid log(0)
    p_nz = p[p > 0]
    return -np.sum(p_nz * np.log2(p_nz))

def run_simulation():
    """
    Simulates the evolution of the state distribution and its entropy
    over multiple policy iterations.
    """
    # --- Simulation Parameters ---
    NUM_STATES = 10
    NUM_ITERATIONS = 25
    # The learning rate alpha models how aggressively the new policy
    # pursues the novel states.
    ALPHA = 0.5
    # Small epsilon to prevent division by zero or log(0).
    EPSILON = 1e-9

    print("This simulation demonstrates how the policy evolves over iterations.\n")

    # --- Initialization (k=0) ---
    # Start with a highly skewed state distribution for the initial policy pi^0.
    # This represents a policy that is initially stuck in one part of the state space.
    p = np.zeros(NUM_STATES)
    p[0] = 1.0
    
    # Calculate and print the initial entropy for pi^0
    entropy_k = calculate_entropy(p)
    max_entropy = np.log2(NUM_STATES)
    print(f"Theoretical Maximum Entropy for {NUM_STATES} states is {max_entropy:.4f}\n")

    print(f"--- Iteration k=0 (Policy pi^0) ---")
    print(f"The initial state distribution is highly skewed.")
    print(f"Entropy H(s) for pi^0: {entropy_k:.4f}\n")

    # --- Iterative Process ---
    # We track the evolution of the state distribution p_k as k increases.
    for k in range(1, NUM_ITERATIONS + 1):
        # The reward for policy pi^k is r_k(s) = -log(p_{k-1}(s)).
        # The new policy pi^k will favor states with high reward. A simple model
        # for the resulting state distribution p_k(s) is one that moves towards
        # a distribution inversely proportional to the previous one.
        
        # This is the target distribution, proportional to exp(reward) = exp(-log(p)) = 1/p
        p_target = 1.0 / (p + EPSILON)
        p_target_normalized = p_target / np.sum(p_target)
        
        # The new distribution is a mix of the old one and the target one.
        p = (1 - ALPHA) * p + ALPHA * p_target_normalized

        # Calculate entropy of the new distribution p_k induced by policy pi^k
        entropy_k = calculate_entropy(p)

        # Print the progress at specific iterations
        if k in [1, 2, 5, 10, 20, 25]:
            print(f"--- Iteration k={k} (Policy pi^{k}) ---")
            print(f"The policy explores, making the state distribution more uniform.")
            print(f"Entropy H(s) for pi^{k}: {entropy_k:.4f}\n")

    print("--- Conclusion ---")
    print(f"As k becomes large, the entropy approaches its maximum value of {max_entropy:.4f}.")
    print("This corresponds to a uniform state distribution.")
    print("Therefore, the policy that maximizes the state entropy is lim_{k->inf} pi^k.")

if __name__ == '__main__':
    run_simulation()