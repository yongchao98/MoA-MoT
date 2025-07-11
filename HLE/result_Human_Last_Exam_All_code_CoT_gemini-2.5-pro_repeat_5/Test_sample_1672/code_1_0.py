import numpy as np

def run_entropy_maximization_simulation():
    """
    Simulates the iterative policy optimization process and tracks state entropy
    to demonstrate which policy maximizes it.
    """
    # --- 1. Environment and RL Parameters ---
    NUM_STATES = 5
    GAMMA = 0.95        # Discount factor
    BETA = 2.0          # Softmax temperature for policy update
    ITERATIONS = 15     # Number of policy iterations k
    VI_ITERATIONS = 100 # Value Iteration steps
    EPSILON = 1e-9      # Small constant to avoid log(0)

    # --- 2. Helper Functions ---
    def get_next_state(s, action_idx):
        """Deterministic transitions for the 1D chain environment."""
        if action_idx == 0:  # Action 'L'
            return max(0, s - 1)
        if action_idx == 1:  # Action 'R'
            return min(NUM_STATES - 1, s + 1)

    def get_transition_matrix(policy):
        """Computes the state transition matrix for a given policy."""
        M = np.zeros((NUM_STATES, NUM_STATES))
        for s in range(NUM_STATES):
            prob_right = policy[s]
            prob_left = 1.0 - prob_right
            M[s, get_next_state(s, 1)] += prob_right
            M[s, get_next_state(s, 0)] += prob_left
        return M

    def get_stationary_distribution(transition_matrix):
        """Computes the stationary distribution of a Markov chain."""
        eigenvalues, eigenvectors = np.linalg.eig(transition_matrix.T)
        closest_eig_idx = np.argmin(np.abs(eigenvalues - 1.0))
        stationary_vector = np.real(eigenvectors[:, closest_eig_idx])
        return stationary_vector / np.sum(stationary_vector)

    def get_entropy(dist):
        """Computes the entropy of a probability distribution in bits."""
        return -np.sum(dist * np.log2(dist + EPSILON))

    # --- 3. Main Simulation Loop ---
    # Start with a biased initial policy pi^0 (heavily favors going right)
    # policy[s] = P(action=R | s)
    pi_k = np.full(NUM_STATES, 0.99)

    max_entropy = np.log2(NUM_STATES)
    print("This simulation shows that the policy converges towards maximum state entropy.")
    print(f"Environment: {NUM_STATES}-state 1D chain.")
    print(f"Maximum possible entropy for this environment: log2({NUM_STATES}) = {max_entropy:.4f}\n")
    print(f"{'Iteration k':<12} | {'State Entropy H(s)':<22} | {'State Distribution p(s)'}")
    print("-" * 80)

    for k in range(ITERATIONS):
        # a. Get state distribution and entropy for the current policy pi_k
        M_k = get_transition_matrix(pi_k)
        p_k = get_stationary_distribution(M_k)
        H_k = get_entropy(p_k)

        # The state distribution p(s) provides the numbers for the entropy equation H(s) = -sum(p*log(p))
        dist_str = np.array2string(p_k, precision=3, floatmode='fixed', suppress_small=True)
        print(f"k = {k:<9} | H(s) = {H_k:<18.4f} | p(s) = {dist_str}")

        # b. Update policy for the next iteration, pi_{k+1}
        # 1. Define intrinsic reward r_{k+1}(s) = -log(p_k(s))
        rewards = -np.log(p_k + EPSILON)

        # 2. Run Value Iteration to find the optimal value function V* for these rewards
        V = np.zeros(NUM_STATES)
        for _ in range(VI_ITERATIONS):
            v_next_L = V[[get_next_state(s, 0) for s in range(NUM_STATES)]]
            v_next_R = V[[get_next_state(s, 1) for s in range(NUM_STATES)]]
            q_values_L = rewards + GAMMA * v_next_L
            q_values_R = rewards + GAMMA * v_next_R
            V = np.maximum(q_values_L, q_values_R)

        # 3. Recalculate Q-values based on the final V*
        v_next_L = V[[get_next_state(s, 0) for s in range(NUM_STATES)]]
        v_next_R = V[[get_next_state(s, 1) for s in range(NUM_STATES)]]
        q_L = rewards + GAMMA * v_next_L
        q_R = rewards + GAMMA * v_next_R
        
        # 4. Update policy pi_{k+1} using a softmax over the Q-values
        exp_q_L = np.exp(BETA * q_L)
        exp_q_R = np.exp(BETA * q_R)
        pi_k = exp_q_R / (exp_q_L + exp_q_R)

if __name__ == '__main__':
    run_entropy_maximization_simulation()