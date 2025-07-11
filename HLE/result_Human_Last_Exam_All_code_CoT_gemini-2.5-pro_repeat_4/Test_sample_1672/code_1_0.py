import numpy as np
import matplotlib.pyplot as plt

def calculate_entropy(p):
    """Calculates the Shannon entropy of a probability distribution."""
    # Filter out zero probabilities to avoid log(0)
    p_nz = p[p > 0]
    return -np.sum(p_nz * np.log2(p_nz))

def simulate_policy_iteration(initial_p, num_iterations):
    """
    Simulates the iterative process of updating the state distribution.
    
    Args:
        initial_p (np.array): The initial state probability distribution.
        num_iterations (int): The number of iterations to run.
    """
    num_states = len(initial_p)
    current_p = np.copy(initial_p)
    
    print("This simulation demonstrates the convergence of the state distribution to a uniform (maximum entropy) distribution.")
    print("-" * 80)
    
    # The target is the uniform distribution
    uniform_dist = np.full(num_states, 1.0 / num_states)
    max_entropy = calculate_entropy(uniform_dist)
    print(f"Target uniform distribution: {np.round(uniform_dist, 3)}")
    print(f"Maximum possible entropy: {max_entropy:.4f}\n")

    print(f"{'Iteration':<10} | {'Current Distribution':<50} | {'Entropy':<10}")
    print("-" * 80)

    for k in range(num_iterations + 1):
        # Print current state
        entropy = calculate_entropy(current_p)
        print(f"{k:<10} | {str(np.round(current_p, 3)):<50} | {entropy:<10.4f}")
        
        if k == num_iterations:
            break

        # The reward is -log(p_k-1(s)). The best response is to go to the state
        # with the minimum probability p_k-1(s).
        min_prob_state_index = np.argmin(current_p)
        
        # The best-response distribution is a delta function on the least likely state.
        best_response_p = np.zeros(num_states)
        best_response_p[min_prob_state_index] = 1.0
        
        # We use a learning rate to ensure convergence, a common practice.
        # A simple choice is 1/(k+1)
        learning_rate = 1.0 / (k + 2) 
        
        # Update the distribution
        current_p = (1 - learning_rate) * current_p + learning_rate * best_response_p

# --- Simulation Parameters ---
# Define a non-uniform initial distribution over 5 states
initial_distribution = np.array([0.5, 0.25, 0.15, 0.05, 0.05])

# Number of iterations to simulate
iterations = 20

# Run the simulation
simulate_policy_iteration(initial_distribution, iterations)