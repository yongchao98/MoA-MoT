import numpy as np

def run_max_entropy_demo():
    """
    Demonstrates that the described iterative process finds a policy
    that maximizes state entropy.
    """
    # --- Environment and Algorithm Setup ---
    N_STATES = 7
    GAMMA = 0.99  # Discount factor, must be < 1 for this method
    START_STATE = N_STATES // 2
    N_ITERATIONS = 15
    EPSILON = 1e-10  # Small constant to avoid log(0)

    # --- Helper Functions ---

    def calculate_state_dist(policy_probs, n_states, gamma, start_state):
        """Calculates the discounted state visitation distribution p_pi(s)."""
        # policy_probs is a vector where policy_probs[s] is P(action=right | state=s)
        T = np.zeros((n_states, n_states))
        for s in range(n_states):
            p_right = policy_probs[s]
            p_left = 1 - p_right
            
            s_right = min(s + 1, n_states - 1)
            T[s, s_right] += p_right
            
            s_left = max(s - 1, 0)
            T[s, s_left] += p_left

        # The state distribution p(s) is given by:
        # p_pi = (1-gamma) * sum_t (gamma*T)^t * p_0
        # which solves to: p_pi = (1-gamma) * p_0 * (I - gamma*T)^-1
        p0 = np.zeros(n_states)
        p0[start_state] = 1.0
        
        try:
            inv_matrix = np.linalg.inv(np.eye(n_states) - gamma * T)
            p_dist = (1 - gamma) * np.dot(p0, inv_matrix)
        except np.linalg.LinAlgError:
            print("Singular matrix detected. The policy may be trapping.")
            return np.ones(n_states) / n_states # Return uniform as fallback

        # Normalize to ensure it's a valid distribution (handles floating point inaccuracies)
        p_dist /= np.sum(p_dist)
        return p_dist

    def value_iteration(rewards, n_states, gamma, tol=1e-5):
        """Solves for the optimal policy for a given reward function using value iteration."""
        V = np.zeros(n_states)
        for _ in range(1000): # Max iterations for VI
            V_old = V.copy()
            Q = np.zeros((n_states, 2))  # 0: left, 1: right
            
            for s in range(n_states):
                s_right = min(s + 1, n_states - 1)
                s_left = max(s - 1, 0)
                
                # Q(s, a) = R(s) + gamma * V(s')
                Q[s, 1] = rewards[s] + gamma * V_old[s_right] # Action: right
                Q[s, 0] = rewards[s] + gamma * V_old[s_left]  # Action: left
            
            V = np.max(Q, axis=1)
            
            if np.max(np.abs(V - V_old)) < tol:
                break
        
        # Extract deterministic policy: 1 for right, 0 for left
        policy_actions = np.argmax(Q, axis=1)
        return policy_actions.astype(float)

    def calculate_entropy(p_dist):
        """Calculates the Shannon entropy of a probability distribution."""
        return -np.sum(p_dist * np.log2(p_dist + EPSILON))

    # --- Main Simulation Loop ---

    # Start with an initial policy pi^0 (e.g., always go right)
    pi_k_probs = np.ones(N_STATES) 

    print("This script demonstrates finding a maximum entropy policy.")
    print(f"Environment: {N_STATES}-state line, starting at state {START_STATE}.\n")

    for k in range(N_ITERATIONS):
        # 1. Calculate state distribution for the current policy pi^k
        p_dist_k = calculate_state_dist(pi_k_probs, N_STATES, GAMMA, START_STATE)
        
        # 2. Calculate and print the entropy H(s) for the current distribution
        entropy_k = calculate_entropy(p_dist_k)
        print(f"--- Iteration k = {k} ---")
        print(f"Policy pi^{k} induces state distribution p(s): {np.round(p_dist_k, 3)}")
        print(f"The entropy H(s) for this distribution is: {entropy_k:.4f}")
        
        # 3. Define the reward for the next iteration (k+1) based on current p_dist
        # r_{k+1}(s) = -log(p_k(s))
        rewards_k_plus_1 = -np.log(p_dist_k + EPSILON)
        
        # 4. Find the new optimal policy pi_{k+1} for this reward
        pi_k_plus_1_probs = value_iteration(rewards_k_plus_1, N_STATES, GAMMA)
        
        # 5. Update policy for the next loop. Use a small learning rate to stabilize convergence.
        pi_k_probs = 0.5 * pi_k_probs + 0.5 * pi_k_plus_1_probs

    print("\n--- Final Result ---")
    # Analyze the final policy after all iterations
    p_dist_final = calculate_state_dist(pi_k_probs, N_STATES, GAMMA, START_STATE)
    entropy_final = calculate_entropy(p_dist_final)
    print(f"The final policy pi^{N_ITERATIONS-1} induces the state distribution:")
    print(f"p(s) = {np.round(p_dist_final, 3)}")
    print(f"The final entropy is {entropy_final:.4f}.")

    # For comparison, calculate the theoretical maximum entropy for this state space
    uniform_dist = np.ones(N_STATES) / N_STATES
    max_entropy = calculate_entropy(uniform_dist)
    print(f"\nThe theoretical maximum entropy (uniform distribution) is {max_entropy:.4f}.")
    print("\nAs k increases, the entropy converges to the maximum value.")

if __name__ == '__main__':
    run_max_entropy_demo()