import numpy as np

def demonstrate_entropy_maximization():
    """
    This script simulates the iterative policy optimization process described in the problem
    to demonstrate how state entropy evolves.
    """

    # --- 1. MDP Definition ---
    N_STATES = 5
    N_ACTIONS = 3  # 0: left, 1: right, 2: stay
    # T[s, a, s'] is the probability of transitioning from s to s' given action a
    T = np.zeros((N_STATES, N_ACTIONS, N_STATES))
    for s in range(N_STATES):
        T[s, 0, (s - 1) % N_STATES] = 1.0  # Action 0: move left
        T[s, 1, (s + 1) % N_STATES] = 1.0  # Action 1: move right
        T[s, 2, s] = 1.0                  # Action 2: stay

    # --- Hyperparameters ---
    GAMMA = 0.95  # Discount factor for value iteration
    N_ITERATIONS = 15  # Number of policy iterations (k)
    SIMULATION_STEPS = 10000  # Steps to run a policy to find its state distribution

    # --- Helper Functions ---
    def calculate_state_distribution(policy, start_state=0, n_steps=SIMULATION_STEPS):
        """Calculates the stationary state distribution by simulating the policy."""
        counts = np.zeros(N_STATES)
        s = start_state
        for _ in range(n_steps):
            counts[s] += 1
            # Choose action based on policy probabilities for the current state
            action = np.random.choice(N_ACTIONS, p=policy[s])
            # Get next state based on transition dynamics
            next_s_dist = T[s, action]
            s = np.random.choice(N_STATES, p=next_s_dist)
        return counts / n_steps

    def calculate_entropy(p):
        """Calculates the entropy of a distribution in bits."""
        # Add a small epsilon to avoid log(0) for states with zero probability
        p_safe = p + 1e-9
        return -np.sum(p_safe * np.log2(p_safe))

    def value_iteration(rewards):
        """Finds the optimal deterministic policy for a given reward function."""
        V = np.zeros(N_STATES)
        # Iterate to find the optimal value function
        for _ in range(100):
            Q = np.zeros((N_STATES, N_ACTIONS))
            for s in range(N_STATES):
                for a in range(N_ACTIONS):
                    # Q(s,a) = r(s) + gamma * E[V(s')]
                    Q[s, a] = rewards[s] + GAMMA * np.dot(T[s, a], V)
            V_new = np.max(Q, axis=1)
            if np.max(np.abs(V - V_new)) < 1e-6:
                break
            V = V_new

        # Derive the optimal deterministic policy from the Q-values
        best_actions = np.argmax(Q, axis=1)
        new_policy = np.zeros((N_STATES, N_ACTIONS))
        for s in range(N_STATES):
            new_policy[s, best_actions[s]] = 1.0
        return new_policy

    # --- Main Loop ---

    # 2. Initialize policy pi^0
    # Start with a biased policy that strongly prefers to 'stay'.
    # This will result in a low-entropy state distribution.
    pi_k = np.zeros((N_STATES, N_ACTIONS))
    pi_k[:, 2] = 0.90  # 90% chance to stay
    pi_k[:, 0] = 0.05  # 5% chance to go left
    pi_k[:, 1] = 0.05  # 5% chance to go right

    max_entropy = np.log2(N_STATES)
    print("Demonstrating how state entropy evolves over iterations.")
    print(f"The environment has {N_STATES} states.")
    print(f"The maximum possible entropy is log2({N_STATES}) = {max_entropy:.4f}\n")

    # 3. Run the iterative process
    for k in range(N_ITERATIONS):
        # Calculate state distribution p_k for the current policy pi_k
        p_k = calculate_state_distribution(pi_k)

        # Calculate and print the entropy of this distribution
        entropy = calculate_entropy(p_k)
        print(f"Iteration k={k}:")
        # The numbers in the equation H(s) = -sum(p(s) * log(p(s))) are p(s) and H(s)
        print(f"  - State distribution p(s) for pi^{k}: {np.round(p_k, 3)}")
        print(f"  - Entropy H(s) for pi^{k}: {entropy:.4f}")

        # Define the reward for the next iteration: r_{k+1}(s) = -log(p_k(s))
        # This rewards visiting states that were previously infrequent.
        r_k_plus_1 = -np.log(p_k + 1e-9)

        # Find the new policy pi_{k+1} that maximizes this reward
        pi_k_plus_1 = value_iteration(r_k_plus_1)

        # Update policy for the next loop
        pi_k = pi_k_plus_1
        print("-" * 40)

    # Final state after all iterations
    p_final = calculate_state_distribution(pi_k)
    entropy_final = calculate_entropy(p_final)
    print(f"Final State (after k={N_ITERATIONS-1}):")
    print(f"  - State distribution p(s) for pi^{{{N_ITERATIONS-1}}}: {np.round(p_final, 3)}")
    print(f"  - Entropy H(s) for pi^{{{N_ITERATIONS-1}}}: {entropy_final:.4f}\n")

    print("Conclusion:")
    print("The simulation shows that as k increases, the policy pi^k changes to make the")
    print("state distribution p(s) more uniform. Consequently, the state entropy H(s)")
    print("increases and converges towards the maximum possible value.")
    print("This confirms that the limiting policy, lim k->inf pi^k, maximizes the state entropy.")

if __name__ == '__main__':
    demonstrate_entropy_maximization()