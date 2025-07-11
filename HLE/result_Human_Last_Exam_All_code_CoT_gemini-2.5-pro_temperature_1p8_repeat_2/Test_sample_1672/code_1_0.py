import numpy as np

def calculate_entropy(p):
    """Calculates the entropy of a probability distribution."""
    p = p[p > 0] # Avoid log(0)
    return -np.sum(p * np.log2(p))

def get_state_transition_matrix(n_states, policy):
    """Computes the state transition matrix P(s'|s) under a given policy."""
    T = np.zeros((n_states, n_states))
    # Actions: 0 for LEFT, 1 for RIGHT
    for s in range(n_states):
        # Action LEFT
        s_next_left = max(0, s - 1)
        T[s, s_next_left] += policy[s, 0]
        # Action RIGHT
        s_next_right = min(n_states - 1, s + 1)
        T[s, s_next_right] += policy[s, 1]
    return T

def get_stationary_distribution(transition_matrix):
    """Calculates the stationary distribution of a Markov chain."""
    eigenvalues, eigenvectors = np.linalg.eig(transition_matrix.T)
    # The stationary distribution is the eigenvector corresponding to the eigenvalue 1.
    stationary_vector = np.real(eigenvectors[:, np.isclose(eigenvalues, 1)])
    # Normalize to get a probability distribution
    return stationary_vector[:, 0] / np.sum(stationary_vector[:, 0])

def value_iteration(n_states, rewards, gamma, n_iters=100):
    """Performs value iteration to find the optimal value function."""
    V = np.zeros(n_states)
    for _ in range(n_iters):
        V_new = np.zeros(n_states)
        for s in range(n_states):
            # Q-values for actions LEFT and RIGHT
            q_left = rewards[s] + gamma * V[max(0, s - 1)]
            q_right = rewards[s] + gamma * V[min(n_states - 1, s + 1)]
            V_new[s] = max(q_left, q_right)
        V = V_new
    return V

def extract_policy(n_states, V, rewards, gamma, temperature=0.1):
    """Extracts a stochastic policy from the value function using softmax."""
    policy = np.zeros((n_states, 2))
    for s in range(n_states):
        # Calculate Q-values from the optimal value function
        q_left = rewards[s] + gamma * V[max(0, s - 1)]
        q_right = rewards[s] + gamma * V[min(n_states - 1, s + 1)]
        
        # Softmax to get a stochastic policy
        q_values = np.array([q_left, q_right]) / temperature
        exp_q = np.exp(q_values - np.max(q_values)) # for numerical stability
        policy[s, :] = exp_q / np.sum(exp_q)
    return policy

def main():
    # --- Parameters ---
    N_STATES = 10
    GAMMA = 0.9
    N_ITERATIONS = 15
    LOG_EPSILON = 1e-10 # To prevent log(0)

    print("Simulating policy iteration with intrinsic motivation reward.")
    print(f"Number of states = {N_STATES}, Max Entropy = log2({N_STATES}) = {np.log2(N_STATES):.4f}\n")

    # --- Initialization (k=0) ---
    # Start with a biased policy: always go right
    policy_k = np.zeros((N_STATES, 2))
    policy_k[:, 1] = 1.0  # Always 100% chance to go RIGHT

    for k in range(N_ITERATIONS):
        # 1. Get state distribution for current policy pi^(k)
        T_k = get_state_transition_matrix(N_STATES, policy_k)
        p_k = get_stationary_distribution(T_k)

        # Calculate and print the entropy
        entropy = calculate_entropy(p_k)
        print(f"Iteration k={k}:")
        print(f"  Policy      = pi^{k}")
        print(f"  Entropy H(s) = {entropy:.4f}")

        if k < N_ITERATIONS - 1:
            # 2. Define reward for next iteration: r_{k+1}(s) = -log(p_k(s))
            # The problem uses r_k = -log(p_{k-1}), so at iteration k we use p_{k-1}
            # Our loop index matches the policy index, so this is correct.
            rewards_k_plus_1 = -np.log(p_k + LOG_EPSILON)

            # 3. Find new policy pi^(k+1) that maximizes this reward
            V_k_plus_1 = value_iteration(N_STATES, rewards_k_plus_1, GAMMA)
            policy_k_plus_1 = extract_policy(N_STATES, V_k_plus_1, rewards_k_plus_1, GAMMA)
            
            # Update policy for the next loop
            policy_k = policy_k_plus_1
        print("-" * 30)
    
    print("As k increases, the entropy H(s) approaches the maximum possible value,")
    print("demonstrating that the state distribution becomes uniform.")
    print("This confirms that the policy that maximizes the entropy is lim_{k->inf} pi^k.")


if __name__ == '__main__':
    main()
