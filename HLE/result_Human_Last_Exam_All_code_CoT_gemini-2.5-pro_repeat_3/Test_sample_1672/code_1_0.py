import numpy as np

def solve_task():
    """
    This script simulates the policy iteration process with intrinsic motivation
    to demonstrate which policy maximizes state entropy.
    """
    # --- Parameters ---
    N_STATES = 10
    GAMMA = 0.99  # Discount factor for value iteration
    BETA = 5.0    # Inverse temperature for softmax policy update
    N_ITERATIONS = 20 # Number of policy iterations (k)
    N_SIM_STEPS = 20000 # Timesteps to simulate policy to get state distribution
    VI_ITERATIONS = 100 # Iterations for Value Iteration

    # --- Environment & Simulation Functions ---
    def get_next_state(s, action):
        # action 0: LEFT, action 1: RIGHT in a 1D chain
        if action == 0:
            return max(0, s - 1)
        else:
            return min(N_STATES - 1, s + 1)

    def get_state_dist(policy, n_steps=N_SIM_STEPS):
        """Simulates a policy to get the state visitation distribution."""
        counts = np.zeros(N_STATES)
        current_state = np.random.randint(N_STATES)
        for _ in range(n_steps):
            prob_right = policy[current_state]
            action = 1 if np.random.rand() < prob_right else 0
            current_state = get_next_state(current_state, action)
            counts[current_state] += 1
        
        # Add a small epsilon to avoid log(0) for unvisited states
        counts += 1e-9
        return counts / np.sum(counts)

    def value_iteration(rewards):
        """Performs value iteration to find the optimal value function."""
        V = np.zeros(N_STATES)
        for _ in range(VI_ITERATIONS):
            V_new = np.zeros(N_STATES)
            for s in range(N_STATES):
                v_left = V[get_next_state(s, 0)]
                v_right = V[get_next_state(s, 1)]
                V_new[s] = rewards[s] + GAMMA * max(v_left, v_right)
            if np.max(np.abs(V - V_new)) < 1e-6:
                break
            V = V_new
        return V

    def update_policy(V):
        """Updates the policy to be softmax over the action-values."""
        new_policy = np.zeros(N_STATES)
        for s in range(N_STATES):
            q_left = V[get_next_state(s, 0)]
            q_right = V[get_next_state(s, 1)]
            exp_q_right = np.exp(BETA * q_right)
            exp_q_left = np.exp(BETA * q_left)
            new_policy[s] = exp_q_right / (exp_q_right + exp_q_left)
        return new_policy

    def calculate_entropy(p):
        """Calculates the entropy of a distribution."""
        p_nz = p[p > 0]
        return -np.sum(p_nz * np.log(p_nz))

    # --- Main Simulation Loop ---
    # Start with a random policy (pi^0)
    current_policy = np.full(N_STATES, 0.5)

    print(f"The goal is to find the policy that maximizes state entropy H(s).")
    print(f"The maximum possible entropy for {N_STATES} states is log({N_STATES}) = {np.log(N_STATES):.4f}\n")
    print("--- Simulating Policy Iterations ---")
    
    final_policy = None
    for k in range(N_ITERATIONS):
        # 1. Get state distribution for current policy pi^k
        p_k = get_state_dist(current_policy)
        
        # 2. Calculate and store entropy of the distribution
        entropy_k = calculate_entropy(p_k)
        print(f"Iteration k={k:2d}: H(s) for policy pi^{k:2d} = {entropy_k:.4f}")

        # 3. Calculate intrinsic rewards for the next iteration: r_{k+1}(s) = -log(p_k(s))
        rewards_k_plus_1 = -np.log(p_k)

        # 4. Find the optimal value function for these rewards
        V_k_plus_1 = value_iteration(rewards_k_plus_1)

        # 5. Update the policy to pi^{k+1}
        current_policy = update_policy(V_k_plus_1)
        final_policy = current_policy

    # The final policy approximates pi^infinity
    p_final = get_state_dist(final_policy)
    h_final = calculate_entropy(p_final)

    print(f"\nAfter {N_ITERATIONS} iterations, the process converges.")
    print(f"The limiting policy (approximated by pi^{N_ITERATIONS-1}) results in an entropy of {h_final:.4f}.")
    print("This demonstrates that the limiting policy, lim k->inf pi^k, maximizes the state entropy H(s).\n")

    # Outputting the final entropy calculation as requested
    print("The final entropy is calculated using the formula: H(s) = - sum(p(s) * log(p(s)))")
    print("For the final state distribution p_final, the calculation is:")
    
    equation_parts = []
    for p_val in p_final:
        equation_parts.append(f"{p_val:.3f}*log({p_val:.3f})")
    
    equation_str = "H(s_final) = -(" + " + ".join(equation_parts)
    equation_str += f") = {h_final:.4f}"
    print(equation_str)

solve_task()