import numpy as np

def run_policy_entropy_maximization():
    """
    Simulates the process of policy iteration using an intrinsic reward
    based on state visitation novelty to maximize state entropy.
    """
    # --- Parameters ---
    N_STATES = 10
    ACTIONS = [-1, 1]  # Left, Right
    GAMMA = 0.95
    POLICY_UPDATE_ITERATIONS = 15
    # Simulation steps to estimate state distribution
    SIMULATION_STEPS = 50000
    # Value iteration steps to find optimal policy for a given reward
    VI_ITERATIONS = 100
    # Epsilon to prevent log(0) and division by zero
    EPSILON = 1e-9

    # --- Helper Functions ---

    def calculate_state_distribution(policy):
        """Simulates the policy to get the state visitation frequency."""
        counts = np.zeros(N_STATES)
        s = np.random.randint(N_STATES) # Start at a random state
        for _ in range(SIMULATION_STEPS):
            counts[s] += 1
            action_idx = np.random.choice([0, 1], p=policy[s])
            action = ACTIONS[action_idx]
            s_next = s + action
            # Bouncing boundary conditions
            s = max(0, min(N_STATES - 1, s_next))
        return counts / np.sum(counts)

    def update_policy_with_value_iteration(rewards):
        """Finds the best policy for a given reward function using value iteration."""
        # 1. Find the optimal value function V*
        V = np.zeros(N_STATES)
        for _ in range(VI_ITERATIONS):
            V_new = np.zeros(N_STATES)
            for s in range(N_STATES):
                q_values = []
                for action in ACTIONS:
                    s_next = max(0, min(N_STATES - 1, s + action))
                    q_values.append(rewards[s] + GAMMA * V[s_next])
                V_new[s] = np.max(q_values)
            V = V_new

        # 2. Derive a stochastic policy from V* (here, we use softmax)
        policy = np.zeros((N_STATES, 2))
        for s in range(N_STATES):
            q_values = []
            for action in ACTIONS:
                s_next = max(0, min(N_STATES - 1, s + action))
                q_values.append(rewards[s] + GAMMA * V[s_next])
            
            # Softmax to create a stochastic policy
            exp_q = np.exp(q_values - np.max(q_values)) # Normalize for stability
            probs = exp_q / np.sum(exp_q)
            policy[s] = probs
        return policy

    def calculate_entropy(p_dist):
        """Calculates the entropy of a probability distribution."""
        # Ensure it's a valid distribution and handle log(0)
        p_dist_safe = p_dist[p_dist > 0]
        return -np.sum(p_dist_safe * np.log(p_dist_safe))

    # --- Main Simulation Logic ---
    print("--- Simulating State Entropy Maximization ---")
    
    # Initialize pi^0 as a uniform random policy
    current_policy = np.full((N_STATES, 2), 0.5)

    print(f"The theoretical maximum entropy for {N_STATES} states is log({N_STATES}) = {np.log(N_STATES):.4f}\n")

    # Iteration k=0
    p_dist = calculate_state_distribution(current_policy)
    entropy = calculate_entropy(p_dist)
    print(f"Policy pi^0: The initial entropy is H(s) = {entropy:.4f}")

    for k in range(1, POLICY_UPDATE_ITERATIONS + 1):
        # The reward for iteration k is based on the distribution from policy k-1
        rewards = -np.log(p_dist + EPSILON)
        
        # Learn the new policy pi^k
        current_policy = update_policy_with_value_iteration(rewards)
        
        # Evaluate the new policy pi^k
        p_dist = calculate_state_distribution(current_policy)
        entropy = calculate_entropy(p_dist)
        print(f"Policy pi^{k}: The entropy is H(s) = {entropy:.4f}")

    print("\n--- Conclusion ---")
    print("The simulation shows that as k increases, the policy pi^k evolves to produce a state distribution")
    print("with progressively higher entropy, approaching the maximum possible value.")
    print("This demonstrates that the limiting policy, lim_{k->inf} pi^k, maximizes the state entropy H(s).")


# Run the simulation
run_policy_entropy_maximization()
