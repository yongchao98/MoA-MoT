import numpy as np

def simulate():
    """
    Simulates the intrinsic motivation process to show entropy maximization.
    """
    # Environment settings
    N_STATES = 10
    ACTIONS = [-1, 1]  # Left, Right
    
    # RL settings
    GAMMA = 0.95
    N_ITERATIONS = 20  # Number of policy updates (k)
    N_SIM_STEPS = 5000 # Steps to run policy to estimate state distribution
    N_VALUE_ITER = 100 # Iterations for Value Iteration

    def transition(s, a):
        """Deterministic transition function for the grid world."""
        return np.clip(s + a, 0, N_STATES - 1)

    def get_state_distribution(policy, start_state=0):
        """
        Runs the policy for N_SIM_STEPS to get the state visitation distribution.
        """
        counts = np.zeros(N_STATES)
        s = start_state
        for _ in range(N_SIM_STEPS):
            action_idx = np.random.choice(len(ACTIONS), p=policy[s])
            a = ACTIONS[action_idx]
            s = transition(s, a)
            counts[s] += 1
        # Add 1 to all counts to avoid zero probabilities (smoothing)
        counts += 1
        return counts / np.sum(counts)

    def update_policy(rewards):
        """
        Updates the policy using value iteration based on the current rewards.
        """
        # 1. Value Iteration to find the optimal value function V*(s)
        V = np.zeros(N_STATES)
        for _ in range(N_VALUE_ITER):
            V_new = np.zeros(N_STATES)
            for s in range(N_STATES):
                q_values = []
                for a in ACTIONS:
                    s_prime = transition(s, a)
                    # The reward is received upon entering the next state s_prime
                    q_val = rewards[s_prime] + GAMMA * V[s_prime]
                    q_values.append(q_val)
                V_new[s] = np.max(q_values)
            V = V_new

        # 2. Policy Extraction from V*(s)
        new_policy = np.zeros((N_STATES, len(ACTIONS)))
        for s in range(N_STATES):
            q_values = []
            for a in ACTIONS:
                s_prime = transition(s, a)
                q_val = rewards[s_prime] + GAMMA * V[s_prime]
                q_values.append(q_val)
            
            # Softmax to convert Q-values to a stochastic policy
            # This encourages exploration within the RL update itself
            q_values = np.array(q_values)
            exp_q = np.exp(q_values - np.max(q_values))
            probs = exp_q / np.sum(exp_q)
            new_policy[s, :] = probs
        return new_policy

    # Initialize a heavily biased policy pi_0 (always go right)
    policy_k = np.zeros((N_STATES, len(ACTIONS)))
    policy_k[:, 1] = 1.0  # 100% probability of going right

    print(f"Simulating for {N_ITERATIONS} iterations in a world with {N_STATES} states.")
    max_entropy = -np.log(1.0 / N_STATES)
    print(f"Maximum possible entropy: {max_entropy:.4f}\n")

    for k in range(N_ITERATIONS):
        # Get state distribution for the current policy
        p_k = get_state_distribution(policy_k)
        
        # Calculate and print entropy
        entropy = -np.sum(p_k * np.log(p_k))
        print(f"Entropy at iteration k={k}: {entropy:.4f}")
        
        # Define rewards for the next iteration
        # Add a small epsilon for numerical stability
        rewards_k_plus_1 = -np.log(p_k + 1e-9)
        
        # Update policy to pi_{k+1}
        policy_k = update_policy(rewards_k_plus_1)

simulate()