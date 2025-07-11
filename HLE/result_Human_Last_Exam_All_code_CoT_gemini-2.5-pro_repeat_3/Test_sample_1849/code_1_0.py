import numpy as np

def run_value_iteration_demo(rewards, gamma, iterations):
    """
    Demonstrates value iteration convergence on a simple 2-state MDP.

    Args:
        rewards (list): A list of two reward values [R(s0), R(s1)].
        gamma (float): The discount factor.
        iterations (int): The number of iterations to run.
    """
    # V is the value function for the two states, initialized to zeros.
    v = np.array([0.0, 0.0])
    
    # In this simple MDP, there is only one action, and each state transitions to itself.
    # V_k+1(s) = R(s) + gamma * V_k(s)
    
    print(f"--- Running demo with rewards {rewards} and gamma={gamma} ---")
    print("Iter | Max Value Change | Change Ratio (-> gamma)")
    print("-----------------------------------------------------")
    
    last_delta = 0
    for i in range(1, iterations + 1):
        v_old = v.copy()
        
        # Apply Bellman backup for each state
        v[0] = rewards[0] + gamma * v_old[0]
        v[1] = rewards[1] + gamma * v_old[1]
        
        # Calculate the max change (infinity norm of the difference)
        delta = np.max(np.abs(v - v_old))
        
        # The ratio of consecutive changes should approach gamma
        ratio = delta / last_delta if last_delta != 0 else float('nan')
        
        # Print the "equation" showing the change at each step
        print(f"{i:4} | {delta:16.4f} | {ratio:.4f}")
        
        last_delta = delta
        # Stop if converged
        if delta < 1e-6:
            print("Converged.")
            break
    print("\n")


# Main execution
if __name__ == "__main__":
    # The discount factor must be < 1 for convergence
    discount_factor = 0.9
    num_iterations = 15

    # Case 1: Rewards in a "standard" range
    small_rewards = [1, -1]
    run_value_iteration_demo(small_rewards, discount_factor, num_iterations)

    # Case 2: Rewards with a large magnitude
    # The convergence rate should be identical, only the values themselves are scaled.
    large_rewards = [1000, -1000]
    run_value_iteration_demo(large_rewards, discount_factor, num_iterations)
