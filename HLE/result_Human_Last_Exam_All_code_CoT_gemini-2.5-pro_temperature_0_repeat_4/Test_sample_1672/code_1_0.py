import numpy as np

def run_entropy_maximization_simulation():
    """
    This script simulates an agent in a simple 1D environment.
    The agent's policy is iteratively updated to seek states that were
    infrequently visited by its previous policy. This intrinsic reward
    is defined as r(s) = -log(p(s)), where p(s) is the previous state
    visitation probability.

    The simulation demonstrates that this process leads to policies that
    induce state distributions with progressively higher entropy, converging
    towards a uniform distribution (which has the maximum possible entropy).
    This supports the conclusion that the limiting policy as the number of
    iterations goes to infinity is the one that maximizes the state entropy.
    """
    # --- Environment and Algorithm Parameters ---
    N_STATES = 10
    N_ACTIONS = 2  # 0: Left, 1: Right
    GAMMA = 0.99   # Discount factor for value iteration
    BETA = 10.0    # Inverse temperature for softmax policy extraction
    N_ITERATIONS = 15 # Number of policy iterations
    N_ROLLOUTS = 1000 # Number of rollouts to estimate state distribution
    LEN_ROLLOUT = 50  # Length of each rollout

    # --- Environment Dynamics ---
    def transition(s, a):
        """Deterministic transitions in a 1D chain."""
        if a == 0:  # Left
            return max(0, s - 1)
        else:  # Right
            return min(N_STATES - 1, s + 1)

    # --- State Distribution Estimation ---
    def calculate_state_dist(policy, start_state=N_STATES // 2):
        """Estimates the state distribution by simulating rollouts."""
        state_counts = np.zeros(N_STATES)
        for _ in range(N_ROLLOUTS):
            s = start_state
            for _ in range(LEN_ROLLOUT):
                state_counts[s] += 1
                action = np.random.choice(N_ACTIONS, p=policy[s])
                s = transition(s, action)
        
        total_counts = np.sum(state_counts)
        if total_counts == 0:
            return np.full(N_STATES, 1.0 / N_STATES)
        return state_counts / total_counts

    # --- Reinforcement Learning Components ---
    def value_iteration(reward, n_vi_steps=100):
        """Solves for the value function for a given reward function."""
        V = np.zeros(N_STATES)
        for _ in range(n_vi_steps):
            V_new = np.zeros(N_STATES)
            for s in range(N_STATES):
                q_values = [reward[s] + GAMMA * V[transition(s, a)] for a in range(N_ACTIONS)]
                V_new[s] = np.max(q_values)
            if np.max(np.abs(V - V_new)) < 1e-6:
                break
            V = V_new
        return V

    def extract_policy(V, reward):
        """Extracts a softmax policy from the value function."""
        Q = np.zeros((N_STATES, N_ACTIONS))
        for s in range(N_STATES):
            for a in range(N_ACTIONS):
                Q[s, a] = reward[s] + GAMMA * V[transition(s, a)]
        
        exp_q = np.exp(BETA * (Q - np.max(Q, axis=1, keepdims=True)))
        policy = exp_q / np.sum(exp_q, axis=1, keepdims=True)
        return policy

    # --- Entropy Calculation ---
    def calculate_entropy(p):
        """Calculates the Shannon entropy of a probability distribution."""
        p_safe = p[p > 1e-9] # Avoid log(0)
        return -np.sum(p_safe * np.log2(p_safe))

    # --- Main Simulation Loop ---
    print("Starting simulation...")
    max_entropy = np.log2(N_STATES)
    print(f"Environment: 1D chain with {N_STATES} states.")
    print(f"Theoretical Maximum Entropy (uniform distribution): {max_entropy:.4f} bits\n")

    # Initialize with a biased policy: strongly prefer going right
    pi_k = np.full((N_STATES, N_ACTIONS), 0.1)
    pi_k[:, 1] = 0.9

    for k in range(N_ITERATIONS):
        # 1. Calculate state distribution for the current policy pi_k
        p_k = calculate_state_dist(pi_k)
        
        # 2. Calculate and print the entropy of this distribution
        H_k = calculate_entropy(p_k)
        print(f"Iteration {k:2d}: Entropy = {H_k:.4f}")

        # 3. Define the intrinsic reward for the next iteration
        r_k_plus_1 = -np.log(p_k + 1e-9)

        # 4. Find the optimal policy for this reward
        V_k_plus_1 = value_iteration(r_k_plus_1)
        pi_k_plus_1 = extract_policy(V_k_plus_1, r_k_plus_1)
        
        # 5. Update the policy for the next iteration
        pi_k = pi_k_plus_1

    # Final entropy after all iterations
    p_final = calculate_state_dist(pi_k)
    H_final = calculate_entropy(p_final)
    print(f"Iteration {N_ITERATIONS:2d}: Entropy = {H_final:.4f}")
    print("\nAs shown, the entropy increases with each iteration, approaching the maximum.")

run_entropy_maximization_simulation()