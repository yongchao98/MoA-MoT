import numpy as np

def simulate():
    """
    Simulates the iterative policy optimization process to show entropy maximization.
    """
    # --- Environment and Simulation Parameters ---
    n_states = 10       # Number of states in a 1D grid
    actions = [-1, 1]   # -1 for Left, 1 for Right
    gamma = 0.99        # Discount factor for value iteration
    n_policy_iter = 15  # Number of policy iterations (k)
    vi_iter = 100       # Steps for value iteration convergence
    sim_steps = 50000   # Simulation steps to find state distribution
    epsilon = 1e-10     # Small number to avoid log(0)

    # --- Helper Functions ---
    def transition(s, a):
        """Deterministic state transition function."""
        return np.clip(s + a, 0, n_states - 1)

    def get_state_distribution(policy, start_state=0):
        """Estimates state distribution by simulating a long trajectory."""
        counts = np.zeros(n_states)
        s = start_state
        for _ in range(sim_steps):
            # Choose action based on the policy (handling deterministic case)
            if np.random.rand() < policy[s, 0]:
                a = actions[0]
            else:
                a = actions[1]
            s = transition(s, a)
            counts[s] += 1
        return counts / sim_steps

    def calculate_entropy(p_dist):
        """Calculates the entropy of a distribution."""
        # Filter out zero probabilities to avoid log(0)
        p_dist_nz = p_dist[p_dist > 0]
        return -np.sum(p_dist_nz * np.log(p_dist_nz))

    # --- Initialization ---
    # Start with a biased policy pi^0 (always goes right)
    policy = np.zeros((n_states, 2))
    policy[:, 1] = 1.0  # Action 1 (Right) has probability 1.0

    print("Starting simulation of iterative policy update.")
    print(f"The maximum possible entropy for {n_states} states is log({n_states}) = {np.log(n_states):.4f}\n")

    # --- Main Iteration Loop ---
    for k in range(n_policy_iter + 1):
        # 1. Get state distribution for the current policy
        p_dist = get_state_distribution(policy)

        # 2. Calculate and print the entropy of this distribution
        entropy = calculate_entropy(p_dist)
        
        print(f"--- Iteration k={k} ---")
        p_dist_str = ", ".join([f"{p:.3f}" for p in p_dist])
        print(f"State Distribution p(s): [{p_dist_str}]")

        # Create the string for the entropy calculation formula
        entropy_terms = [f"({p:.3f} * log({p:.3f}))" for p in p_dist if p > epsilon]
        entropy_eq_str = " + ".join(entropy_terms)
        print(f"Entropy H(p) = -[ {entropy_eq_str} ] = {entropy:.4f}\n")

        if k == n_policy_iter:
            print("Simulation finished. Final policy leads to the highest entropy.")
            break
        
        # 3. Calculate rewards for the next iteration: r_{k+1}(s) = -log(p_k(s))
        rewards = -np.log(p_dist + epsilon)

        # 4. Find the next policy pi^{k+1} using Value Iteration
        V = np.zeros(n_states)
        for _ in range(vi_iter):
            Q = np.zeros((n_states, len(actions)))
            for s in range(n_states):
                for i, a in enumerate(actions):
                    s_next = transition(s, a)
                    # The reward is associated with arriving in the next state
                    Q[s, i] = rewards[s_next] + gamma * V[s_next]
            V_new = np.max(Q, axis=1)
            if np.max(np.abs(V - V_new)) < 1e-5:
                break
            V = V_new
        
        # 5. Derive the new greedy policy
        Q_final = np.zeros((n_states, len(actions)))
        for s in range(n_states):
            for i, a in enumerate(actions):
                s_next = transition(s, a)
                Q_final[s, i] = rewards[s_next] + gamma * V[s_next]

        best_actions = np.argmax(Q_final, axis=1)
        
        # Update policy to be deterministic based on the best action
        new_policy = np.zeros((n_states, 2))
        for s in range(n_states):
            new_policy[s, best_actions[s]] = 1.0
        
        policy = new_policy

if __name__ == '__main__':
    simulate()
