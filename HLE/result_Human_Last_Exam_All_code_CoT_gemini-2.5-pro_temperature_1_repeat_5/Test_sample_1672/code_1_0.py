import numpy as np

def calculate_entropy(p):
    """Calculates the entropy of a probability distribution using log base 2."""
    # Add a small epsilon to avoid log(0) for states with zero probability
    epsilon = 1e-12
    # Ensure probabilities are clipped between epsilon and 1.0
    p_clipped = np.clip(p, epsilon, 1.0)
    return -np.sum(p_clipped * np.log2(p_clipped))

def run_simulation(num_states=10, num_iterations=15, alpha=0.5):
    """
    Simulates the evolution of the state distribution towards maximum entropy.

    This simulation models the core idea of the problem. At each step k,
    the new state distribution p_k is encouraged to visit states that were
    infrequent under the previous distribution p_{k-1}.

    Args:
        num_states (int): The number of states in the environment.
        num_iterations (int): The number of policy iterations k.
        alpha (float): The learning rate for updating the distribution. A value
                       between 0 and 1 that controls how fast the distribution
                       moves towards the target.
    """
    print(f"Starting simulation with {num_states} states...")
    print("The goal is to see how the state distribution entropy changes over iterations.")
    
    # 1. Initialize p_0(s), the state distribution for the initial policy pi^0.
    # We create a skewed, non-uniform distribution.
    p = np.full(num_states, 1.0)
    p[0] = num_states * 2  # Make the first state much more likely
    p = p / np.sum(p)  # Normalize to make it a valid probability distribution
    
    print("\n--- Iteration k=0 (Initial Policy pi^0) ---")
    print(f"Initial state distribution p_0(s): {np.round(p, 4)}")
    entropy = calculate_entropy(p)
    print(f"Entropy H(s) at k=0: {entropy:.4f}")
    
    # The maximum possible entropy for this number of states (uniform distribution)
    max_entropy = np.log2(num_states)
    print(f"Maximum possible entropy for {num_states} states is: {max_entropy:.4f}\n")

    # 2. Run the iterative process for k=1, 2, ...
    for k in range(1, num_iterations + 1):
        # The reward for the next policy is r_{k}(s) = -log(p_{k-1}(s)).
        # This reward is high for states 's' where p_{k-1}(s) is low.
        # The new policy pi^k will learn to visit these high-reward states.
        # We model the resulting state distribution p_k(s) by moving it towards a
        # target distribution 'q' that is inversely proportional to the current 'p'.
        
        # Define the target distribution q_k(s) ~ 1 / p_{k-1}(s)
        epsilon = 1e-12
        inverse_p = 1 / (p + epsilon)
        q = inverse_p / np.sum(inverse_p)

        # Update the state distribution: p_k = (1-alpha)*p_{k-1} + alpha*q_k
        p = (1 - alpha) * p + alpha * q
        
        print(f"--- Iteration k={k} ---")
        entropy = calculate_entropy(p)
        print(f"Updated state distribution p_{k}(s): {np.round(p, 4)}")
        print(f"Entropy H(s) at k={k}: {entropy:.4f}")
        
    print("\n--- Simulation End ---")
    print(f"After {num_iterations} iterations, the distribution is nearly uniform.")
    print(f"The final entropy ({entropy:.4f}) is very close to the maximum possible entropy ({max_entropy:.4f}).")
    print("This demonstrates that as k increases, the policy pi^k induces a state distribution with increasing entropy,")
    print("converging to the maximum entropy distribution in the limit.")

if __name__ == '__main__':
    run_simulation()