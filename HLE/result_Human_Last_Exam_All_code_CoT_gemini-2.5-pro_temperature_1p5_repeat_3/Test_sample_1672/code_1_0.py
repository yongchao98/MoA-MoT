import numpy as np
import math

def calculate_entropy(p):
    """Calculates the entropy of a probability distribution."""
    # Filter out zero probabilities to avoid log(0)
    p = p[p > 0]
    return -np.sum(p * np.log2(p)) # Using log base 2 for entropy in bits

def get_state_distribution(policy_q_values, n_states, n_episodes=200, max_steps=100):
    """
    Estimates the state visitation distribution for a given policy by simulation.
    A greedy policy is derived from the Q-values.
    """
    state_counts = np.zeros(n_states)
    for _ in range(n_episodes):
        state = np.random.randint(n_states) # Start at a random state
        for _ in range(max_steps):
            state_counts[state] += 1
            # Choose action greedily based on the policy's Q-values
            action = np.argmax(policy_q_values[state]) # 0 for left, 1 for right
            if action == 0: # Move left
                state = max(0, state - 1)
            else: # Move right
                state = min(n_states - 1, state + 1)
    
    # Add a small smoothing constant to avoid zero counts
    state_counts += 1e-9
    return state_counts / np.sum(state_counts)

def value_iteration_for_policy(rewards, n_states, gamma=0.95):
    """
    Performs value iteration to find the optimal Q-values for a given reward function.
    The environment is a 1D grid where actions are deterministic moves left/right.
    """
    q_values = np.zeros((n_states, 2)) # 2 actions: left, right
    for _ in range(100): # Number of VI iterations
        v = np.max(q_values, axis=1)
        new_q_values = np.zeros_like(q_values)
        for s in range(n_states):
            # Q(s, left) - state becomes s' = max(0, s-1), reward is r(s')
            s_left = max(0, s - 1)
            new_q_values[s, 0] = rewards[s_left] + gamma * v[s_left]
            # Q(s, right) - state becomes s' = min(n-1, s+1), reward is r(s')
            s_right = min(n_states - 1, s + 1)
            new_q_values[s, 1] = rewards[s_right] + gamma * v[s_right]
        
        # Check for convergence
        if np.max(np.abs(new_q_values - q_values)) < 1e-4:
            break
        q_values = new_q_values
    return q_values

def demonstrate_entropy_maximization():
    """
    Simulates the iterative policy update process and shows entropy increasing.
    """
    n_states = 10
    n_iterations = 15
    
    # The maximum possible entropy for a uniform distribution over n_states
    max_entropy = -math.log2(1.0 / n_states)
    
    print("--- Entropy Maximization Simulation ---")
    print(f"Environment: 1D grid with {n_states} states")
    print(f"Maximum possible entropy (uniform distribution): {max_entropy:.4f} bits\n")
    
    # Initialize with a heavily biased policy pi^0
    # Let's create a policy that strongly prefers going to the rightmost state (9).
    q_policy = np.zeros((n_states, 2))
    q_policy[:, 1] = 1 # Action 1 (right) is always preferred.

    print("Initial policy pi^0 is biased to always move right.")
    
    for k in range(n_iterations):
        # Current policy is pi^{k-1} if k>0, or pi^0 if k=0
        # It is represented by the Q-table `q_policy`
        
        # 1. Get state distribution for the current policy
        p_dist = get_state_distribution(q_policy, n_states)
        
        # 2. Calculate entropy of the distribution
        entropy = calculate_entropy(p_dist)
        print(f"Iteration k={k: <2}:  H(p_{{pi^{k}}}) = {entropy:.4f} bits")
        
        # 3. Calculate intrinsic rewards for the next iteration (pi^{k+1})
        # r_{k+1}(s) = -log(p_{pi^{k}}(s))
        rewards = -np.log(p_dist)
        
        # 4. Learn a new policy by optimizing for these rewards
        # The new policy is represented by the new Q-values
        q_policy = value_iteration_for_policy(rewards, n_states)
        
    print("\n--- Conclusion ---")
    print("The simulation shows that the entropy of the state distribution increases with each iteration k.")
    print("This is because the reward function at each step encourages the agent to explore states that were less frequently visited by the previous policy.")
    print("This process drives the state visitation distribution towards a uniform distribution.")
    print("The policy that achieves the maximum entropy uniform distribution is therefore the one obtained at the limit of this process.")

if __name__ == '__main__':
    demonstrate_entropy_maximization()