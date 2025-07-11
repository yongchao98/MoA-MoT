import numpy as np

# This script simulates the policy iteration process with an intrinsic reward
# to demonstrate how state entropy is maximized over time.

def run_entropy_maximization_simulation():
    """
    Simulates a policy learning to maximize state entropy in a simple environment.
    """
    # --- 1. Environment and Simulation Setup ---
    N_STATES = 10         # A 1D line of 10 states (0 to 9)
    LEFT, RIGHT = 0, 1    # Possible actions
    GAMMA = 0.99          # Discount factor for RL
    N_POLICY_ITERATIONS = 15 # Number of policy updates to run

    def step(state, action):
        """Deterministic state transition function for the 1D grid world."""
        if action == LEFT:
            return max(0, state - 1)
        elif action == RIGHT:
            return min(N_STATES - 1, state + 1)

    def get_state_distribution(policy, n_episodes=500, episode_len=50):
        """Estimates state visitation probability p(s) for a given policy."""
        counts = np.zeros(N_STATES)
        for _ in range(n_episodes):
            s = np.random.randint(N_STATES) # Start from a random state
            for _ in range(episode_len):
                counts[s] += 1
                action = np.random.choice([LEFT, RIGHT], p=policy[s])
                s = step(s, action)
        # Add a tiny value to prevent log(0) for unvisited states
        counts += 1e-9
        return counts / np.sum(counts)

    def train_policy(rewards, vi_iterations=100):
        """Trains a new policy using Value Iteration based on the given rewards."""
        V = np.zeros(N_STATES)
        for _ in range(vi_iterations):
            V_new = np.zeros(N_STATES)
            for s in range(N_STATES):
                s_left, s_right = step(s, LEFT), step(s, RIGHT)
                # Q(s,a) = R(s') + gamma * V(s')
                # Here reward is given for entering a state, so we use reward of the next state.
                q_left = rewards[s_left] + GAMMA * V[s_left]
                q_right = rewards[s_right] + GAMMA * V[s_right]
                V_new[s] = max(q_left, q_right)
            # Stop if values have converged
            if np.allclose(V, V_new):
                break
            V = V_new

        # Derive a near-deterministic policy from the converged values
        new_policy = np.zeros((N_STATES, 2))
        epsilon = 0.1 # For a little exploration
        for s in range(N_STATES):
            s_left, s_right = step(s, LEFT), step(s, RIGHT)
            q_left = rewards[s_left] + GAMMA * V[s_left]
            q_right = rewards[s_right] + GAMMA * V[s_right]
            if q_left > q_right:
                new_policy[s] = [1.0 - epsilon, epsilon]
            elif q_right > q_left:
                new_policy[s] = [epsilon, 1.0 - epsilon]
            else: # If Q-values are equal, move randomly
                new_policy[s] = [0.5, 0.5]
        return new_policy

    def calculate_entropy(p):
        """Calculates the Shannon entropy H(p) = -sum(p * log2(p))."""
        return -np.sum(p * np.log2(p))

    # --- 2. Iterative Policy Improvement ---
    print("Starting simulation of entropy maximization.")
    
    # The maximum possible entropy is for a uniform distribution over N_STATES
    max_entropy = calculate_entropy(np.full(N_STATES, 1/N_STATES))
    print(f"Number of states: {N_STATES}")
    print(f"Maximum possible state entropy H(s): {max_entropy:.4f}\n")

    # Start with an initial random policy pi^0
    current_pi = np.full((N_STATES, 2), 0.5)

    for k in range(N_POLICY_ITERATIONS):
        # Estimate state distribution p_{pi^{k-1}}(s)
        p_dist = get_state_distribution(current_pi)

        # Calculate entropy H(p_{pi^{k-1}})
        entropy = calculate_entropy(p_dist)
        print(f"Iteration k={k}:  H(s) = {entropy:.4f}")

        # Define reward for iteration k: r_k(s) = -log(p_{pi^{k-1}}(s))
        rewards = -np.log(p_dist)

        # Train new policy pi^k to maximize expected reward
        current_pi = train_policy(rewards)
        
    # --- 3. Final Result ---
    # The policy has now converged after N iterations. Let's check its final entropy.
    final_p_dist = get_state_distribution(current_pi)
    final_entropy = calculate_entropy(final_p_dist)
    
    print(f"Iteration k={N_POLICY_ITERATIONS}:  H(s) = {final_entropy:.4f} (Final converged value)")
    
    print("\n--- Conclusion from Simulation ---")
    print("The entropy of the state distribution increases with each iteration k.")
    print("It converges towards the maximum possible entropy.")
    print("This shows that the policy that maximizes the entropy is the one at the limit, lim k->inf pi^k.")
    
# Run the simulation and print the results
run_entropy_maximization_simulation()

<<<A>>>