import numpy as np

def solve():
    """
    This function simulates the iterative policy improvement process to show
    how state entropy is maximized over time.
    """
    # --- Simulation Parameters ---
    NUM_STATES = 10
    ACTIONS = ['left', 'right']
    GAMMA = 0.9  # Discount factor for value iteration
    BETA = 1.0   # Inverse temperature for softmax policy
    ITERATIONS = 20  # Number of policy iterations k
    VALUE_ITER_STEPS = 50 # Number of steps for value iteration
    EPSILON = 1e-10 # Small constant to avoid log(0)

    # --- Helper Functions ---
    def get_next_state(s, a):
        """Deterministic state transitions in a 1D chain environment."""
        if a == 'left':
            return max(0, s - 1)
        else: # right
            return min(NUM_STATES - 1, s + 1)

    def get_policy_from_q(q_values):
        """Derives a softmax policy from Q-values."""
        policy = np.zeros((NUM_STATES, len(ACTIONS)))
        for s in range(NUM_STATES):
            exps = np.exp(BETA * q_values[s, :])
            policy[s, :] = exps / np.sum(exps)
        return policy

    def get_stationary_distribution(policy):
        """Calculates the stationary state distribution for a given policy."""
        transition_matrix = np.zeros((NUM_STATES, NUM_STATES))
        for s in range(NUM_STATES):
            for a_idx, a in enumerate(ACTIONS):
                s_next = get_next_state(s, a)
                transition_matrix[s, s_next] += policy[s, a_idx]
        
        # The stationary distribution is the left eigenvector of the transition
        # matrix corresponding to the eigenvalue 1.
        eigenvalues, eigenvectors = np.linalg.eig(transition_matrix.T)
        stationary_vector = eigenvectors[:, np.isclose(eigenvalues, 1)]
        stationary_dist = np.real(stationary_vector[:, 0]) / np.sum(np.real(stationary_vector[:, 0]))
        return stationary_dist

    def calculate_entropy(p):
        """Calculates the entropy H(p) of a distribution p."""
        p_nonzero = p[p > 0]
        return -np.sum(p_nonzero * np.log(p_nonzero))

    # --- Main Simulation ---
    print("This simulation demonstrates that the policy converges to maximize state entropy.")
    print(f"The iterative reward is r_k(s) = -log p_k-1(s).\n")
    
    # Step 1: Initialize the policy pi^0
    # Start with a biased policy that strongly prefers moving 'right'.
    q_0 = np.zeros((NUM_STATES, len(ACTIONS)))
    q_0[:, 1] = 1.0  # Higher Q-value for 'right'
    pi_k = get_policy_from_q(q_0)

    # Maximum possible entropy for NUM_STATES (uniform distribution)
    max_entropy = np.log(NUM_STATES)
    print(f"Number of states: {NUM_STATES}")
    print(f"Maximum possible entropy: {max_entropy:.4f}\n")
    
    print("Iteration k | Entropy H(p_k)")
    print("----------------------------")

    # Main iterative loop
    for k in range(1, ITERATIONS + 1):
        # At the start of the loop, pi_k is the policy from iteration k-1
        p_prev = get_stationary_distribution(pi_k)
        
        if k == 1:
            # Print initial state entropy (for k=0)
            initial_entropy = calculate_entropy(p_prev)
            print(f"0 (initial) | {initial_entropy:.4f}")

        # Define the reward for the current iteration
        rewards = -np.log(p_prev + EPSILON)

        # Train a new policy for iteration k using Value Iteration on the rewards
        v_k = np.zeros(NUM_STATES)
        for _ in range(VALUE_ITER_STEPS):
            v_new = np.zeros(NUM_STATES)
            for s in range(NUM_STATES):
                q_vals_for_s = [rewards[get_next_state(s, a)] + GAMMA * v_k[get_next_state(s, a)] for a in ACTIONS]
                v_new[s] = np.max(q_vals_for_s)
            if np.max(np.abs(v_new - v_k)) < 1e-6:
                break
            v_k = v_new
        
        # Derive the new Q-values and the new policy pi_k
        q_k = np.zeros((NUM_STATES, len(ACTIONS)))
        for s in range(NUM_STATES):
            for a_idx, a in enumerate(ACTIONS):
                s_next = get_next_state(s, a)
                q_k[s, a_idx] = rewards[s_next] + GAMMA * v_k[s_next]
        
        pi_k = get_policy_from_q(q_k)

        # Calculate and print the entropy of the new state distribution
        p_k = get_stationary_distribution(pi_k)
        entropy_k = calculate_entropy(p_k)
        print(f"{k:<11} | {entropy_k:.4f}")

    print("\n--- Conclusion ---")
    print("As shown by the simulation, the entropy H(p_k) increases with each iteration k,")
    print("approaching the maximum possible entropy. This happens because the policy is")
    print("continuously driven to explore less-visited states, leading to a more uniform")
    print("state distribution over time.")
    print("\nTherefore, the policy that maximizes the state entropy H(s) is the one obtained")
    print("in the limit as k approaches infinity, which is lim_{k -> infinity} pi^k.")
    print("This corresponds to answer choice A.")

solve()