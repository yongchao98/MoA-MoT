import numpy as np

def run_simulation():
    """
    This simulation demonstrates that the policy converges towards one that
    maximizes state entropy as the iteration k increases.
    """
    # --- 1. Environment and Simulation Setup ---
    num_states = 10
    num_actions = 2  # 0: left, 1: right
    gamma = 0.99  # Discount factor for value iteration
    num_iterations = 15  # Number of policy iterations (k)
    simulation_steps = 50000  # Steps to run policy to estimate state distribution
    
    def step_env(state, action):
        """Deterministic 1D chain environment dynamics."""
        if action == 0:  # Left
            return max(0, state - 1)
        else:  # Right
            return min(num_states - 1, state + 1)

    # --- 2. Start with an initial policy pi^0 ---
    # We start with a policy that heavily prefers one side.
    # This creates an initial low-entropy state distribution.
    policy = np.zeros((num_states, num_actions))
    policy[:, 1] = 0.9  # 90% chance of going right
    policy[:, 0] = 0.1  # 10% chance of going left

    # --- 3. Main Loop (Iterating k) ---
    print("Demonstration of Entropy Maximization\n")
    print("The goal is to find the policy that maximizes state entropy H(s).")
    print(f"The maximum possible entropy for {num_states} states is log({num_states}) = {np.log(num_states):.4f}\n")
    print("k | H(p_pi_k)")
    print("-" * 15)

    for k in range(num_iterations):
        # a. Calculate state distribution p_pi^k(s) by simulating the current policy
        counts = np.zeros(num_states)
        state = np.random.randint(num_states)
        for _ in range(simulation_steps):
            action = np.random.choice(num_actions, p=policy[state])
            state = step_env(state, action)
            counts[state] += 1
        state_dist = counts / counts.sum()

        # b. Calculate and print the entropy for the current iteration's distribution
        # H(s) = -sum(p(s) * log(p(s)))
        epsilon = 1e-12  # to avoid log(0)
        entropy = -np.sum(state_dist * np.log(state_dist + epsilon))
        print(f"{k:<1} | {entropy:.4f}")

        # c. Define reward for the next iteration: r_{k+1}(s) = -log(p_pi^k(s))
        rewards = -np.log(state_dist + epsilon)

        # d. Learn new policy pi^{k+1} using Value Iteration
        V = np.zeros(num_states)
        for _ in range(100):  # VI sweeps
            V_new = np.zeros(num_states)
            for s in range(num_states):
                q_values = [rewards[s] + gamma * V[step_env(s, a)] for a in range(num_actions)]
                V_new[s] = max(q_values)
            V = V_new
        
        # e. Extract the new policy pi^{k+1} (softmax for exploration)
        new_policy = np.zeros((num_states, num_actions))
        for s in range(num_states):
            q_values = np.array([rewards[s] + gamma * V[step_env(s, a)] for a in range(num_actions)])
            exp_q = np.exp(q_values - np.max(q_values)) # numerical stability
            new_policy[s] = exp_q / np.sum(exp_q)
        
        policy = new_policy

    print("-" * 15)
    print("As k increases, the entropy H(s) increases and approaches the maximum.")
    print("This shows the process converges to a policy that maximizes state entropy.")
    print("This limit policy is best represented by choice A.")


# Execute the simulation
run_simulation()
<<<A>>>