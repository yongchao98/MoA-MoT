import numpy as np

def demonstrate_entropy_maximization():
    """
    This simulation demonstrates how an iterative policy update scheme based on
    intrinsic rewards leads to a policy that maximizes state entropy.
    """

    # 1. Environment Setup
    # A simple 1D grid world with 5 states: 0-1-2-3-4
    # Actions: 0 (left), 1 (right)
    N_STATES = 5
    N_ACTIONS = 2
    GAMMA = 0.99 # Discount factor for value iteration

    def get_transitions(s, a):
        """Deterministic transitions: moving left from 0 stays at 0, right from 4 stays at 4."""
        if a == 0:  # left
            return max(0, s - 1)
        if a == 1:  # right
            return min(N_STATES - 1, s + 1)
        return s

    # 2. Helper Functions
    def get_state_dist(policy):
        """Calculates the stationary state distribution for a given policy."""
        P_pi = np.zeros((N_STATES, N_STATES))
        for s_prime in range(N_STATES):
            for action in range(N_ACTIONS):
                s_next = get_transitions(s_prime, action)
                P_pi[s_prime, s_next] += policy[s_prime, action]
        
        # Solve p * P_pi = p for p, which is the left eigenvector of P_pi for eigenvalue 1.
        eigenvalues, eigenvectors = np.linalg.eig(P_pi.T)
        idx = np.argmin(np.abs(eigenvalues - 1.0))
        p = np.real(eigenvectors[:, idx])
        return p / np.sum(p) # Normalize to a probability distribution

    def value_iteration(rewards):
        """Finds the optimal deterministic policy for a given reward function."""
        V = np.zeros(N_STATES)
        for _ in range(1000):
            V_old = V.copy()
            q_values = np.zeros((N_STATES, N_ACTIONS))
            for s in range(N_STATES):
                for a in range(N_ACTIONS):
                    s_next = get_transitions(s, a)
                    # Reward is given for being in state s
                    q_values[s, a] = rewards[s] + GAMMA * V[s_next]
            V = np.max(q_values, axis=1)
            if np.max(np.abs(V - V_old)) < 1e-6:
                break
        
        # Extract policy
        policy = np.zeros((N_STATES, N_ACTIONS))
        best_actions = np.argmax(q_values, axis=1)
        policy[np.arange(N_STATES), best_actions] = 1.0
        return policy

    def get_entropy(p):
        """Calculates the entropy H(p) = -sum(p_i * log(p_i))."""
        p_safe = p[p > 1e-9] # Avoid log(0)
        return -np.sum(p_safe * np.log2(p_safe)) # Use log base 2 for bits

    # 3. Main Simulation Loop
    print("--- Simulating Entropy Maximization ---")
    print(f"The goal is to find the policy that maximizes state entropy H(s).")
    print(f"The policy pi^k is updated using rewards r_k(s) = -log(p_{{pi^(k-1)}}(s)).\n")

    # Start with a non-uniform policy: always go right
    pi_k = np.zeros((N_STATES, N_ACTIONS))
    pi_k[:, 1] = 1.0

    for k in range(10):
        # Calculate state distribution and entropy for the current policy pi^k
        p_k = get_state_dist(pi_k)
        H_k = get_entropy(p_k)

        print(f"--- Iteration k={k} ---")
        print(f"The state distribution p_{{pi^{k}}}(s) is: {np.round(p_k, 3)}")
        print(f"The entropy H(s) for this distribution is: {H_k:.4f} bits")
        
        # Define rewards for the next iteration (k+1)
        # Add a small epsilon to avoid -log(0) for unvisited states
        rewards = -np.log(p_k + 1e-9)
        
        # Find the new policy pi_{k+1} by optimizing for these rewards
        pi_k_plus_1 = value_iteration(rewards)
        
        # Update policy for the next loop
        pi_k = pi_k_plus_1
        print("-" * 30)

    # Final state
    p_final = get_state_dist(pi_k)
    H_final = get_entropy(p_final)
    print("--- Final State (after k=9) ---")
    print(f"The final state distribution is: {np.round(p_final, 3)}")
    print(f"The final entropy is: {H_final:.4f} bits")

    # Compare with maximum possible entropy
    max_entropy = np.log2(N_STATES)
    print(f"\nThe maximum possible entropy for {N_STATES} states is log2({N_STATES}) = {max_entropy:.4f} bits.")
    print("As k increases, the entropy approaches this maximum value.")
    print("This demonstrates that the limit policy, lim_{k->inf} pi^k, maximizes the state entropy.")

if __name__ == '__main__':
    demonstrate_entropy_maximization()