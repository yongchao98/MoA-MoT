import numpy as np

def calculate_entropy(p):
    """Calculates the entropy of a probability distribution."""
    # Add a small epsilon to prevent log(0)
    epsilon = 1e-9
    return -np.sum(p * np.log2(p + epsilon))

def run_simulation():
    """
    Simulates the evolution of the state distribution and its entropy.
    """
    N_STATES = 10
    N_ITERATIONS = 20

    # Start with a highly skewed initial state distribution p_0(s)
    p = np.array([0.5] + [0.5 / (N_STATES - 1)] * (N_STATES - 1))
    
    print(f"Starting simulation with {N_STATES} states.")
    print(f"Maximum theoretical entropy: log2({N_STATES}) = {np.log2(N_STATES):.4f}\n")

    for k in range(N_ITERATIONS):
        # Calculate the entropy of the current state distribution
        entropy = calculate_entropy(p)
        
        # Print the current state
        print(f"--- Iteration k = {k} ---")
        # The prompt requires printing each number in the final equation.
        # Here we show the distribution p_k(s) that goes into the entropy equation H(s) = -sum(p * log(p)).
        print("State distribution p_k(s):")
        print(np.round(p, 4))
        print(f"Entropy H(s): {entropy:.4f}")
        print("-" * (25 + len(str(k))))

        # If entropy is very close to maximum, we have converged
        if np.isclose(entropy, np.log2(N_STATES)):
            print("\nConvergence reached: State distribution is uniform.")
            break

        # Calculate the intrinsic reward r_k+1(s) = -log(p_k(s))
        # Add epsilon to avoid log(0)
        epsilon = 1e-9
        reward = -np.log(p + epsilon)

        # The new policy will drive the agent towards high-reward states.
        # We model the resulting state distribution as being proportional to the reward.
        # This is a simplification of a full reinforcement learning update.
        p_target = reward / np.sum(reward)
        
        # The new state distribution p_k+1 is a mix of the old and the target distributions,
        # representing a single step of policy improvement.
        p = 0.7 * p + 0.3 * p_target

run_simulation()
