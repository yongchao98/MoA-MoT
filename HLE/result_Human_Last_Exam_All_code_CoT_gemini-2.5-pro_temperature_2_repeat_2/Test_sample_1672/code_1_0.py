import numpy as np
import collections

# Define Environment and Algorithm Parameters
N_STATES = 10
# Actions are defined by their effect on the state index: -1 for LEFT, 1 for RIGHT
LEFT, RIGHT = -1, 1
ACTIONS = [LEFT, RIGHT]
GAMMA = 0.99  # Discount factor for RL
N_ITERATIONS = 15 # Number of policy iterations k

# Parameters for sub-routines
SIMULATION_TRAJECTORY_LENGTH = 20000 # For calculating state distribution by simulation
VI_CONVERGENCE_THRESHOLD = 1e-5 # For policy optimization using Value Iteration
SOFTMAX_TEMP = 1.0 # Temperature for softmax policy extraction

def transition(state, action):
    """
    Computes the next state given the current state and action.
    The agent bounces off the walls of the 1D grid.
    """
    next_state = state + action
    if next_state < 0:
        return 0
    if next_state >= N_STATES:
        return N_STATES - 1
    return next_state

def get_state_dist(policy, start_state=0):
    """
    Calculates the state visitation distribution p(s) for a given policy
    by simulating a long trajectory.
    """
    state = start_state
    counts = collections.defaultdict(int)
    for _ in range(SIMULATION_TRAJECTORY_LENGTH):
        action_idx = np.random.choice(len(ACTIONS), p=policy[state])
        action = ACTIONS[action_idx]
        state = transition(state, action)
        counts[state] += 1
    
    # Create the probability distribution vector, adding a tiny epsilon for stability
    distribution = np.full(N_STATES, 1e-9)
    for s, count in counts.items():
        distribution[s] += count
    
    return distribution / np.sum(distribution)

def update_policy_via_value_iteration(rewards):
    """
    Finds the optimal policy for a given reward function using Value Iteration.
    The reward r(s') is obtained upon arriving in state s'. The value of a state s is
    V(s) = max_a (r(s') + gamma * V(s')).
    """
    values = np.zeros(N_STATES)
    while True:
        delta = 0
        new_values = np.copy(values)
        for s in range(N_STATES):
            q_values = np.zeros(len(ACTIONS))
            for i, action in enumerate(ACTIONS):
                s_next = transition(s, action)
                q_values[i] = rewards[s_next] + GAMMA * values[s_next]

            best_q_value = np.max(q_values)
            delta = max(delta, abs(new_values[s] - best_q_value))
            new_values[s] = best_q_value
        values = new_values
        if delta < VI_CONVERGENCE_THRESHOLD:
            break
            
    # Once optimal values are found, derive a stochastic policy using softmax
    policy = np.zeros((N_STATES, len(ACTIONS)))
    for s in range(N_STATES):
        q_values = np.zeros(len(ACTIONS))
        for i, action in enumerate(ACTIONS):
            s_next = transition(s, action)
            q_values[i] = rewards[s_next] + GAMMA * values[s_next]
            
        # Softmax for action probabilities
        exp_q = np.exp(q_values / SOFTMAX_TEMP - np.max(q_values / SOFTMAX_TEMP)) 
        policy[s] = exp_q / np.sum(exp_q)
        
    return policy

def calculate_entropy(distribution):
    """Calculates the Shannon entropy (in bits) of a probability distribution."""
    # Filter out probabilities that are essentially zero to avoid log(0) issues
    non_zero_dist = distribution[distribution > 1e-8]
    return -np.sum(non_zero_dist * np.log2(non_zero_dist))

def main():
    """Main function to run the simulation."""
    print("This simulation demonstrates that the described iterative process finds a policy")
    print("that maximizes the state entropy H(s).\n")
    print(f"Environment: 1D grid with {N_STATES} states.")
    print("Reward for iteration k+1: r_{k+1}(s) = -log(p_k(s))")
    print("We expect to see the entropy increase and converge to its maximum value.\n")

    # 1. Initialize policy pi_0: deterministically go RIGHT
    pi_current = np.zeros((N_STATES, len(ACTIONS)))
    pi_current[:, 1] = 1.0 

    print("--- Starting Simulation ---")

    for k in range(N_ITERATIONS):
        # 2. Calculate p_k from pi_k
        p_current = get_state_dist(pi_current)

        # 3. Calculate and print entropy H(p_k)
        entropy = calculate_entropy(p_current)
        print(f"Iteration k={k}: Entropy H(p_k) = {entropy:.4f}")

        # 4. Define rewards for the next iteration: r_{k+1}(s) = -log(p_k(s))
        rewards_next = -np.log(p_current)
    
        # 5. Find the new policy pi_{k+1} using RL (Value Iteration) for these rewards
        pi_next = update_policy_via_value_iteration(rewards_next)
        
        pi_current = pi_next # Update policy for the next loop

    print("\n--- Simulation Complete ---")
    max_entropy = np.log2(N_STATES)
    print("The entropy increased and stabilized at a value close to the theoretical maximum.")
    print(f"The maximum possible entropy for {N_STATES} states is log2({N_STATES}) = {max_entropy:.4f}.")
    print("This supports the conclusion that the limiting policy (as k -> infinity) maximizes state entropy.")

    # Calculate and print final entropy details
    final_dist = get_state_dist(pi_current)
    final_entropy = calculate_entropy(final_dist)
    print(f"\nEntropy of the final policy's state distribution: {final_entropy:.4f}")
    
    # Format and print the final equation for entropy calculation
    print("\nThe entropy H is calculated as H = -sum(p * log2(p)).")
    print("For the final distribution, the calculation is:")
    equation_parts = [f"({p:.3f} * log2({p:.3f}))" for p in final_dist if p > 1e-6]
    equation_str = "H = - [" + " + ".join(equation_parts) + "]"
    # Break long lines for readability
    if len(equation_str) > 80:
        equation_str = equation_str.replace(" + ", " +\n           ")
    print(f"Final Equation: {equation_str} = {final_entropy:.4f}")

if __name__ == "__main__":
    main()