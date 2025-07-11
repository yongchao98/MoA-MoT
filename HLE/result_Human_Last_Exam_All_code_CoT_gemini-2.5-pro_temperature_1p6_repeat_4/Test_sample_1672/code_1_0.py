import numpy as np

def solve():
    """
    This script demonstrates the concept that an agent whose reward is based
    on the novelty of states will, over time, develop a policy that
    maximizes state entropy. It simulates an agent on a 1D line.
    """

    # --- Environment Parameters ---
    N_STATES = 15  # The number of states in our 1D world
    LEFT = 0
    RIGHT = 1

    # --- Algorithm Parameters ---
    N_ITERATIONS = 20  # Number of policy updates (k)
    N_SIMULATION_STEPS = 50000  # Steps to run to estimate state distribution
    GAMMA = 0.95  # Discount factor for value iteration
    TEMPERATURE = 0.05  # Temperature for softmax policy extraction
    EPSILON = 1e-10  # Small value to avoid log(0)

    def get_state_dist(policy, start_state=N_STATES // 2):
        """Simulates the policy to get the state visitation distribution."""
        counts = np.zeros(N_STATES)
        current_state = start_state
        for _ in range(N_SIMULATION_STEPS):
            counts[current_state] += 1
            # Choose action based on the policy's probability distribution for the current state
            action = np.random.choice([LEFT, RIGHT], p=policy[current_state])
            if action == LEFT:
                current_state = max(0, current_state - 1)
            else:  # RIGHT
                current_state = min(N_STATES - 1, current_state + 1)
        return counts / np.sum(counts)

    def update_policy(rewards):
        """Finds the optimal policy for a given reward function using value iteration."""
        # 1. Value Iteration to find V(s)
        V = np.zeros(N_STATES)
        for _ in range(150):  # Iterate until convergence
            V_old = V.copy()
            Q = np.zeros((N_STATES, 2))
            for s in range(N_STATES):
                next_s_left = max(0, s - 1)
                Q[s, LEFT] = rewards[s] + GAMMA * V_old[next_s_left]

                next_s_right = min(N_STATES - 1, s + 1)
                Q[s, RIGHT] = rewards[s] + GAMMA * V_old[next_s_right]
            
            V = Q.max(axis=1)
            if np.max(np.abs(V - V_old)) < 1e-6:
                break
        
        # 2. Re-calculate Q-values with the final V(s)
        Q_final = np.zeros((N_STATES, 2))
        for s in range(N_STATES):
            next_s_left = max(0, s - 1)
            Q_final[s, LEFT] = rewards[s] + GAMMA * V[next_s_left]
            next_s_right = min(N_STATES - 1, s + 1)
            Q_final[s, RIGHT] = rewards[s] + GAMMA * V[next_s_right]

        # 3. Softmax policy extraction from Q-values
        # Using a softmax function allows for a stochastic policy.
        # A higher temperature leads to more random actions.
        scaled_Q = Q_final / TEMPERATURE
        # Subtract max for numerical stability to avoid overflow in exp
        max_q = np.max(scaled_Q, axis=1, keepdims=True)
        exp_Q = np.exp(scaled_Q - max_q)
        policy = exp_Q / np.sum(exp_Q, axis=1, keepdims=True)
        return policy

    def calculate_entropy(p_dist):
        """Calculates the entropy of a probability distribution."""
        # Filter out zero probabilities to avoid issues with log(0)
        p_dist_nz = p_dist[p_dist > 0]
        return -np.sum(p_dist_nz * np.log2(p_dist_nz))

    # --- Main Simulation Loop ---
    
    # The maximum possible entropy corresponds to a uniform distribution
    max_entropy = np.log2(N_STATES)
    print("This simulation shows that as k increases, the policy pi^k evolves to")
    print("maximize the state entropy H(s). The entropy of the resulting state")
    print("distribution approaches the theoretical maximum, which corresponds")
    print("to a uniform state visitation distribution.\n")
    print(f"Theoretical Maximum Entropy for {N_STATES} states: {max_entropy:.4f}\n")
    print(f"{'Iteration k':<15} | {'State Entropy H(s)':<25}")
    print("-" * 45)

    # Initialize with a uniform random policy (pi^0)
    policy_k = np.full((N_STATES, 2), 0.5)

    for k in range(N_ITERATIONS):
        # 1. Estimate state distribution p_k(s) for the current policy pi^k
        p_k = get_state_dist(policy_k)

        # 2. Calculate and print the entropy of this distribution
        entropy_k = calculate_entropy(p_k)
        print(f"{k:<15} | {entropy_k:<25.4f}")

        # 3. Define the intrinsic rewards for the next iteration (k+1)
        # States that were rare under pi^k get a high reward
        rewards_k_plus_1 = -np.log(p_k + EPSILON)

        # 4. Find the next policy pi^{k+1} that maximizes this reward
        policy_k = update_policy(rewards_k_plus_1)

solve()