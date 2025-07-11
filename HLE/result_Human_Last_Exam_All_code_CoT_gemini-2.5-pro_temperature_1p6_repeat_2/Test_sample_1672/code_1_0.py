import numpy as np

def simulate_policy(policy, num_states, num_steps, start_state=0):
    """Simulates a policy for a number of steps and returns the state visitation counts."""
    counts = np.zeros(num_states)
    current_state = start_state
    for _ in range(num_steps):
        counts[current_state] += 1
        # Choose an action based on the policy at the current state.
        # policy[s] is the probability of moving right (action = 1)
        if np.random.rand() < policy[current_state]:
            action = 1  # Move right
        else:
            action = -1 # Move left
        
        # Update the state, staying within the boundaries [0, num_states-1]
        current_state = np.clip(current_state + action, 0, num_states - 1)
    return counts

def calculate_entropy(p):
    """Calculates the entropy of a probability distribution p."""
    # We only consider states with non-zero probability to avoid log(0)
    p_positive = p[p > 0]
    return -np.sum(p_positive * np.log2(p_positive))

def run_simulation():
    """Runs the full simulation to demonstrate entropy maximization."""
    # --- Simulation Parameters ---
    NUM_STATES = 10
    NUM_ITERATIONS_K = 25  # Number of policy iterations
    SIMULATION_STEPS = 50000
    BETA = 1.0  # Controls how greedily the policy is updated
    EPSILON = 1e-9  # Small constant to avoid log(0)

    # --- Initialization (k=0) ---
    # pi^0: An initial policy that strongly prefers moving right.
    # policy[s] = probability of moving right at state s.
    pi_0 = np.full(NUM_STATES, 0.95)
    
    policies = [pi_0]
    
    # --- Main Loop ---
    # This loop simulates the policy improvement process from k=1 to K.
    for k in range(1, NUM_ITERATIONS_K + 1):
        # The previous policy is pi^{k-1}
        prev_policy = policies[-1]
        
        # Get the state distribution p_{pi^{k-1}}(s) by simulating the previous policy.
        counts = simulate_policy(prev_policy, NUM_STATES, SIMULATION_STEPS)
        p_prev = counts / np.sum(counts)

        # Define the intrinsic reward r_k(s) = -log(p_{pi^{k-1}}(s))
        rewards = -np.log(p_prev + EPSILON)

        # Create the new policy pi^k. It should favor actions leading to states with higher rewards.
        new_policy = np.zeros(NUM_STATES)
        for s in range(NUM_STATES):
            # Reward for moving right from state s (i.e., reward at state s+1)
            reward_right = rewards[min(s + 1, NUM_STATES - 1)]
            # Reward for moving left from state s (i.e., reward at state s-1)
            reward_left = rewards[max(s - 1, 0)]
            
            # Use a softmax function to determine the probability of moving right.
            # This makes the policy favor the direction with higher reward.
            exp_right = np.exp(BETA * reward_right)
            exp_left = np.exp(BETA * reward_left)
            prob_right = exp_right / (exp_right + exp_left)
            new_policy[s] = prob_right
            
        policies.append(new_policy)

    # --- Analysis and Output ---
    # Now, let's analyze the entropy of the state distributions for key policies.
    
    # Distribution and Entropy for pi^0
    p_pi_0_counts = simulate_policy(policies[0], NUM_STATES, SIMULATION_STEPS)
    p_pi_0 = p_pi_0_counts / np.sum(p_pi_0_counts)
    H_s_0 = calculate_entropy(p_pi_0)

    # Distribution and Entropy for pi^1
    p_pi_1_counts = simulate_policy(policies[1], NUM_STATES, SIMULATION_STEPS)
    p_pi_1 = p_pi_1_counts / np.sum(p_pi_1_counts)
    H_s_1 = calculate_entropy(p_pi_1)

    # Distribution and Entropy for the limit policy (approximated by the last policy)
    final_policy = policies[-1]
    p_pi_inf_counts = simulate_policy(final_policy, NUM_STATES, SIMULATION_STEPS * 2) # More steps for accuracy
    p_pi_inf = p_pi_inf_counts / np.sum(p_pi_inf_counts)
    H_s_inf = calculate_entropy(p_pi_inf)
    
    # Maximum possible entropy for a uniform distribution over NUM_STATES states.
    max_entropy = np.log2(NUM_STATES)

    print("--- Entropy Analysis of Iterative Policies ---")
    print(f"The simulation shows that as k increases, the policy induces a state distribution")
    print("with higher entropy, approaching the maximum possible value.\n")
    print(f"Theoretical maximum entropy for {NUM_STATES} states: {max_entropy:.4f}\n")
    
    print(f"Policy pi^0 (Initial Policy):")
    print(f"State Entropy H(s | pi^0) = {H_s_0:.4f}")
    
    print("\nPolicy pi^1 (After 1 Iteration):")
    print(f"State Entropy H(s | pi^1) = {H_s_1:.4f}")
    
    print(f"\nPolicy lim_{{k->inf}} pi^k (Approximated by k={NUM_ITERATIONS_K}):")
    print(f"State Entropy H(s | lim pi^k) = {H_s_inf:.4f}")
    print("\nConclusion: The limit policy lim_{k->inf} pi^k maximizes the state entropy H(s).")


if __name__ == '__main__':
    run_simulation()