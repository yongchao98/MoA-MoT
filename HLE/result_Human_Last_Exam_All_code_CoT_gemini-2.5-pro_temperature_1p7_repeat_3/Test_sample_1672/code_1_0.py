import numpy as np

def run_entropy_maximization_demo():
    """
    This function demonstrates that an iterative policy update, where the reward is
    the negative log-probability of the previous policy's state distribution,
    leads to a policy that maximizes state entropy.
    """
    N_STATES = 10
    
    # --- Helper Functions ---
    
    def get_transition_matrix(policy):
        """
        Calculates the environment's state transition matrix given a policy.
        - policy[s]: probability of moving 'right' from state s.
        - 1-policy[s]: probability of moving 'left'.
        The world is a 1D grid that wraps around (a cycle).
        """
        T = np.zeros((N_STATES, N_STATES))
        for s in range(N_STATES):
            p_right = policy[s]
            p_left = 1 - p_right
            s_right = (s + 1) % N_STATES
            s_left = (s - 1 + N_STATES) % N_STATES
            T[s, s_right] += p_right
            T[s, s_left] += p_left
        return T

    def get_stationary_distribution(T):
        """
        Calculates the stationary state distribution p(s) for a transition matrix T.
        This corresponds to the eigenvector of T.T with eigenvalue 1.
        """
        eigenvalues, eigenvectors = np.linalg.eig(T.T)
        idx = np.argmin(np.abs(eigenvalues - 1.0))
        stationary_dist = np.real(eigenvectors[:, idx])
        # Normalize to ensure it's a valid probability distribution
        return stationary_dist / np.sum(stationary_dist)

    def calculate_entropy(p):
        """Calculates the entropy of a distribution: H(p) = -sum(p * log2(p))."""
        p_safe = p + 1e-12 # Add a small epsilon to avoid log(0)
        return -np.sum(p_safe * np.log2(p_safe))

    def update_policy(reward, temperature=0.1):
        """
        Updates the policy to prefer actions leading to higher-reward states.
        Uses a softmax function for probabilistic choice.
        """
        new_policy = np.zeros(N_STATES)
        for s in range(N_STATES):
            s_right = (s + 1) % N_STATES
            s_left = (s - 1 + N_STATES) % N_STATES
            
            reward_right = reward[s_right]
            reward_left = reward[s_left]
            
            exp_r_right = np.exp(reward_right / temperature)
            exp_r_left = np.exp(reward_left / temperature)
            
            new_policy[s] = exp_r_right / (exp_r_right + exp_r_left)
        return new_policy

    # --- Main Simulation ---

    # Start with an initial policy pi^0 that is strongly biased.
    # It heavily prefers the 'right' action in all states.
    policy_k = np.full(N_STATES, 0.95)
    
    # The maximum possible entropy is for a uniform distribution: log2(N)
    max_entropy = np.log2(N_STATES)
    
    print("--- Policy Iteration for Entropy Maximization ---")
    print(f"Environment: A {N_STATES}-state cyclic world.")
    print(f"Goal: Maximize state entropy H(s). Theoretical Maximum: {max_entropy:.4f}")
    print("-" * 51)
    
    # Run for 15 iterations
    for k in range(15):
        # For the current policy pi^k, find its induced state distribution p_k
        T_k = get_transition_matrix(policy_k)
        p_k = get_stationary_distribution(T_k)
        
        # Calculate the entropy of this distribution
        h_k = calculate_entropy(p_k)
        
        print(f"Iteration k = {k:2d}: Entropy H(p_pi^{k}(s)) = {h_k:.4f}")

        # The reward for the *next* iteration (k+1) is based on p_k
        reward_k_plus_1 = -np.log(p_k + 1e-9)
        
        # Update the policy from pi^k to pi^{k+1} using this reward
        policy_k = update_policy(reward_k_plus_1)

    print("-" * 51)
    print("Observation: As k increases, the policy is refined to make the state")
    print("distribution more uniform, causing the entropy to approach its maximum.")
    print("This implies the entropy-maximizing policy is the limit as k -> infinity.")


if __name__ == '__main__':
    run_entropy_maximization_demo()