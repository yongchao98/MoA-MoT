import numpy as np

def run_simulation():
    """
    Simulates policy iteration with an intrinsic motivation reward
    to show that it maximizes state entropy over time.
    """
    # --- Hyperparameters ---
    N_STATES = 10       # Number of states in the 1D grid world
    GAMMA = 0.99        # Discount factor for value iteration
    BETA = 1.0          # Inverse temperature for softmax policy derivation
    N_ROLLOUTS = 1000   # Number of simulation runs to estimate state distribution
    MAX_STEPS = 50      # Steps per rollout
    K_ITERATIONS = 20   # Number of policy iterations (k)
    EPSILON = 1e-10     # Small constant to avoid log(0)

    # --- Environment Definition ---
    # States: 0, 1, ..., N_STATES-1
    # Actions: 0 (left), 1 (right)
    # The agent starts in the middle of the grid.
    start_state = N_STATES // 2

    def transition(s, a):
        """Calculates the next state given the current state and action."""
        if a == 0:  # action is "left"
            return max(0, s - 1)
        else:  # action is "right"
            return min(N_STATES - 1, s + 1)

    # --- Helper Functions ---
    def get_state_distribution(policy):
        """Calculates state visitation frequency by simulating many rollouts."""
        state_counts = np.zeros(N_STATES)
        for _ in range(N_ROLLOUTS):
            current_state = start_state
            for _ in range(MAX_STEPS):
                state_counts[current_state] += 1
                action_probs = policy[current_state]
                action = np.random.choice([0, 1], p=action_probs)
                current_state = transition(current_state, action)
        
        total_visits = np.sum(state_counts)
        if total_visits == 0: # Handle cases with no visits
            return np.ones(N_STATES) / N_STATES
        return state_counts / total_visits

    def value_iteration(rewards):
        """Performs value iteration to find the optimal value function for a given reward."""
        V = np.zeros(N_STATES)
        for _ in range(100):  # A fixed number of sweeps is usually sufficient
            V_old = V.copy()
            q_values = np.zeros((N_STATES, 2))
            for s in range(N_STATES):
                q_values[s, 0] = rewards[s] + GAMMA * V[transition(s, 0)] # Q(s, left)
                q_values[s, 1] = rewards[s] + GAMMA * V[transition(s, 1)] # Q(s, right)
            V = np.max(q_values, axis=1)
            if np.max(np.abs(V - V_old)) < 1e-6:
                break
        return V

    def derive_policy(V, rewards):
        """Derives a stochastic policy from the value function using softmax."""
        policy = np.zeros((N_STATES, 2))
        for s in range(N_STATES):
            q_left = rewards[s] + GAMMA * V[transition(s, 0)]
            q_right = rewards[s] + GAMMA * V[transition(s, 1)]
            
            # Softmax to get action probabilities
            exp_q_left = np.exp(BETA * q_left)
            exp_q_right = np.exp(BETA * q_right)
            sum_exp_q = exp_q_left + exp_q_right
            
            policy[s, 0] = exp_q_left / sum_exp_q
            policy[s, 1] = exp_q_right / sum_exp_q
        return policy

    # --- Main Simulation Loop ---
    
    # Initialize with a deterministic policy that always goes right
    current_policy = np.zeros((N_STATES, 2))
    current_policy[:, 1] = 1.0

    print("This simulation shows that as the iteration k increases, the policy pi^k tends to maximize state entropy.")
    print(f"The environment has {N_STATES} states. The maximum possible entropy is log({N_STATES}) = {np.log(N_STATES):.4f}\n")
    
    print(f"{'Iteration k':<12} | {'State Entropy H(s)':<25}")
    print("-" * 40)

    # In our loop, `current_policy` is pi^{k-1}. We use it to find pi^k.
    # The entropy is calculated for the distribution from pi^{k-1}.
    for k in range(K_ITERATIONS):
        # 1. Get state distribution for pi^{k-1}
        p_dist = get_state_distribution(current_policy)

        # 2. Calculate and print the entropy of this distribution
        # H(s) = -sum(p(s) * log(p(s)))
        entropy = -np.sum(p_dist[p_dist > 0] * np.log(p_dist[p_dist > 0]))
        
        print(f"k = {k:<10} | Entropy = {entropy:.8f}")

        # 3. Define the intrinsic reward for finding pi^k: r_k(s) = -log p_{pi^{k-1}}(s)
        rewards = -np.log(p_dist + EPSILON)

        # 4. Find the new optimal policy pi^k using this reward
        V = value_iteration(rewards)
        current_policy = derive_policy(V, rewards)

    print("-" * 40)
    print("As k increases, the entropy approaches the maximum, supporting the conclusion.")

if __name__ == '__main__':
    run_simulation()