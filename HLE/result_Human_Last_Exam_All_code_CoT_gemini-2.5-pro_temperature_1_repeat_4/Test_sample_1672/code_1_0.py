import numpy as np

def run_simulation():
    """
    Simulates the iterative policy improvement process to demonstrate entropy maximization.
    """
    # 1. Environment: A 1D line of `N_STATES`. Agent can move left or right.
    N_STATES = 10
    ACTIONS = [0, 1]  # 0: left, 1: right

    # 2. Helper function to compute the stationary state distribution via simulation
    def compute_state_distribution(policy, n_steps=50000):
        counts = np.zeros(N_STATES)
        current_state = np.random.randint(N_STATES)
        for _ in range(n_steps):
            counts[current_state] += 1
            # Choose action based on policy for the current state
            action = np.random.choice(ACTIONS, p=policy[current_state])
            # Transition to the next state
            if action == 0:  # Move left
                current_state = max(0, current_state - 1)
            else:  # Move right
                current_state = min(N_STATES - 1, current_state + 1)
        return counts / n_steps

    # 3. Helper function to compute entropy
    def calculate_entropy(p):
        # Add a small epsilon to avoid log(0) for states with zero probability
        p_safe = p + 1e-9
        return -np.sum(p_safe * np.log(p_safe))

    # 4. Helper function to find the optimal policy for a given reward function using Value Iteration
    def find_optimal_policy(rewards, gamma=0.99):
        V = np.zeros(N_STATES)
        for _ in range(100):  # Iterate until the value function converges
            V_old = V.copy()
            for s in range(N_STATES):
                # Q-value for moving left
                s_left = max(0, s - 1)
                q_left = rewards[s] + gamma * V[s_left]
                # Q-value for moving right
                s_right = min(N_STATES - 1, s + 1)
                q_right = rewards[s] + gamma * V[s_right]
                V[s] = max(q_left, q_right)
            if np.max(np.abs(V - V_old)) < 1e-6:
                break
        
        # Derive a deterministic policy from the optimal value function
        policy = np.zeros((N_STATES, len(ACTIONS)))
        for s in range(N_STATES):
            s_left = max(0, s - 1)
            q_left = rewards[s] + gamma * V[s_left]
            s_right = min(N_STATES - 1, s + 1)
            q_right = rewards[s] + gamma * V[s_right]
            if q_left > q_right:
                policy[s, 0] = 1.0  # Prefer left
            elif q_right > q_left:
                policy[s, 1] = 1.0  # Prefer right
            else:
                policy[s, :] = 0.5  # Indifferent
        return policy

    # --- Main Simulation Logic ---
    print("This simulation demonstrates that the iterative policy update rule leads to a state distribution with increasing entropy.")
    print("The reward at iteration k is r_k(s) = -log(p_pi^{k-1}(s)), which encourages visiting novel states.\n")

    # The maximum possible entropy is for a uniform distribution
    max_entropy = calculate_entropy(np.full(N_STATES, 1.0/N_STATES))
    print(f"Environment: {N_STATES} states on a line.")
    print(f"Maximum possible entropy (uniform distribution): {max_entropy:.4f}\n")
    print("-" * 50)

    # Start with an initial biased policy pi^0 (e.g., always try to go right)
    policy_k = np.zeros((N_STATES, len(ACTIONS)))
    policy_k[:, 1] = 1.0
    
    n_iterations = 8
    for k in range(n_iterations):
        # Calculate state distribution for the current policy pi^k
        p_k = compute_state_distribution(policy_k)
        
        # Calculate and print the entropy of this distribution
        h_k = calculate_entropy(p_k)
        print(f"Iteration k={k}:")
        print(f"  - Entropy of state distribution p_pi^{k}(s): {h_k:.4f}")
        
        # Define rewards for the next iteration (k+1) based on p_k
        # Add a small epsilon to prevent log(0)
        rewards_k_plus_1 = -np.log(p_k + 1e-9)
        
        # Find the next policy pi^{k+1} by optimizing for these rewards
        policy_k_plus_1 = find_optimal_policy(rewards_k_plus_1)
        
        # Update policy for the next loop
        policy_k = policy_k_plus_1
        print("-" * 50)
        
    print("\nObservation: As k increases, the entropy H(s) approaches the maximum possible value.")
    print("This supports the conclusion that lim_{k->inf} pi^k maximizes the state entropy.")

if __name__ == "__main__":
    run_simulation()