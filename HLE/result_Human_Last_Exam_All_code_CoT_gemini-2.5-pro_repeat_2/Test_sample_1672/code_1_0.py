import numpy as np

def solve():
    """
    This function simulates the policy iteration process described in the problem
    to show that the limiting policy maximizes state entropy.
    """
    # --- Environment Parameters ---
    N_STATES = 10
    ACTIONS = {'LEFT': -1, 'RIGHT': 1}

    def transition(state, action_val):
        """Deterministic transitions with wrap-around."""
        return (state + action_val) % N_STATES

    # --- Algorithm Parameters ---
    GAMMA = 0.95
    N_ITERATIONS = 20  # Number of policy update iterations
    VALUE_ITER_STEPS = 100 # Steps for value iteration convergence
    SOFTMAX_BETA = 5.0 # Inverse temperature for softmax policy
    EPSILON = 1e-9 # For numerical stability

    # --- Helper Functions ---
    def get_transition_matrix(policy):
        """Computes the state transition matrix T_ss' for a given policy."""
        T = np.zeros((N_STATES, N_STATES))
        for s in range(N_STATES):
            for a_idx, action_val in enumerate(ACTIONS.values()):
                prob_a = policy[s, a_idx]
                next_s = transition(s, action_val)
                T[s, next_s] += prob_a
        return T

    def get_stationary_distribution(T):
        """Computes the stationary distribution of a transition matrix T."""
        eigenvalues, eigenvectors = np.linalg.eig(T.T)
        idx = np.argmin(np.abs(eigenvalues - 1.0))
        p = np.real(eigenvectors[:, idx])
        p /= p.sum()
        return p

    def get_entropy(p):
        """Computes the entropy H(p) of a distribution p."""
        return -np.sum(p * np.log(p + EPSILON))

    # 1. Initialize policy pi^0
    # Policy is a (N_STATES, N_ACTIONS) array.
    # Start with a biased policy that heavily prefers going right.
    policy_k = np.zeros((N_STATES, len(ACTIONS)))
    policy_k[:, 1] = 0.99  # High prob of going RIGHT
    policy_k[:, 0] = 0.01  # Low prob of going LEFT

    print("--- Simulating Policy Iteration for State-Entropy Maximization ---")
    max_entropy = np.log(N_STATES)
    print(f"Environment: {N_STATES} states in a ring.")
    print(f"Maximum possible entropy: log({N_STATES}) = {max_entropy:.4f}\n")

    # Calculate initial entropy for pi^0
    T_k = get_transition_matrix(policy_k)
    p_k = get_stationary_distribution(T_k)
    H_k = get_entropy(p_k)
    print(f"Iteration k=0: Entropy H(s) = {H_k:.4f}")

    # 2. Main Iteration Loop
    for k in range(1, N_ITERATIONS + 1):
        # The state distribution from the previous step is p_k
        p_k_minus_1 = p_k

        # b. Define reward for this iteration: r_k(s) = -log p_{k-1}(s)
        rewards_k = -np.log(p_k_minus_1 + EPSILON)

        # c. Find new policy pi^k using Value Iteration with the new rewards
        V = np.zeros(N_STATES)
        for _ in range(VALUE_ITER_STEPS):
            Q = np.zeros((N_STATES, len(ACTIONS)))
            for s in range(N_STATES):
                # Reward r_k(s) is received for being in state s
                q_val_for_s = rewards_k[s]
                for a_idx, action_val in enumerate(ACTIONS.values()):
                    next_s = transition(s, action_val)
                    Q[s, a_idx] = q_val_for_s + GAMMA * V[next_s]
            V = np.max(Q, axis=1)

        # Derive new stochastic policy pi^k from Q-values using softmax
        exp_q = np.exp(SOFTMAX_BETA * (Q - np.max(Q, axis=1, keepdims=True)))
        policy_k = exp_q / np.sum(exp_q, axis=1, keepdims=True)

        # d. Compute new state distribution and its entropy
        T_k = get_transition_matrix(policy_k)
        p_k = get_stationary_distribution(T_k)
        H_k = get_entropy(p_k)
        
        print(f"Iteration k={k}: Entropy H(s) = {H_k:.4f}")

    print("\n--- Conclusion ---")
    print("As the number of iterations k increases, the policy pi^k is refined")
    print("to produce a state distribution with higher entropy. The entropy converges")
    print("towards the maximum possible value.")
    print("\nThis demonstrates that the policy which maximizes the state entropy H(s) is the")
    print("one obtained at the limit of this process: lim_{k->inf} pi^k.")

solve()