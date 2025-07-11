import numpy as np

def demonstrate_value_iteration(reward_value):
    """
    Demonstrates value iteration on a simple 2-state MDP.

    Args:
        reward_value (float): The reward for transitioning from state 0 to state 1.
    """
    # MDP parameters
    num_states = 2
    num_actions = 2
    gamma = 0.9  # Discount factor
    threshold = 1e-5  # Convergence threshold

    # V[s] is the value of state s
    V = np.zeros(num_states)

    # Transition probabilities: P[s, a, s']
    P = np.zeros((num_states, num_actions, num_states))
    # From s0, action 0 -> s1; action 1 -> s0 (stay)
    P[0, 0, 1] = 1.0
    P[0, 1, 0] = 1.0
    # From s1, both actions lead to s1 (stay)
    P[1, 0, 1] = 1.0
    P[1, 1, 1] = 1.0

    # Rewards: R[s, a, s']
    R = np.zeros((num_states, num_actions, num_states))
    # The only non-zero reward is for (s0, a0) -> s1
    R[0, 0, 1] = reward_value

    print(f"--- Running Value Iteration with R(s0,a0->s1) = {reward_value} ---")
    iteration = 0
    while True:
        iteration += 1
        delta = 0
        V_old = V.copy()
        
        # The update equation: V_{k+1}(s) = max_a Σ_{s'} P(s'|s,a)[R(s,a,s') + γV_k(s')]
        for s in range(num_states):
            action_values = []
            for a in range(num_actions):
                # Calculate the value for taking action 'a' in state 's'
                # This is the inner sum: Σ_{s'} P(...)[R(...) + γV(...)]
                expected_value = sum([P[s, a, s_next] * (R[s, a, s_next] + gamma * V_old[s_next])
                                     for s_next in range(num_states)])
                action_values.append(expected_value)
            
            # Update V(s) with the max over action values
            V[s] = max(action_values)

        # Check for convergence
        delta = np.max(np.abs(V - V_old))
        if delta < threshold:
            break
            
    print(f"Convergence reached after {iteration} iterations.")
    print(f"Final Value Function: V(s0) = {V[0]:.2f}, V(s1) = {V[1]:.2f}\n")
    # For this specific MDP, the exact solution is:
    # V*(s1) = 0
    # V*(s0) = max(reward_value + gamma * V*(s1), 0 + gamma * V*(s0))
    # -> V*(s0) = max(reward_value, gamma * V*(s0))
    # This implies V*(s0) = reward_value if reward_value > 0, and V*(s0) = 0 if reward_value <= 0.
    
# Case 1: Large positive reward
demonstrate_value_iteration(1000)

# Case 2: Large negative reward
demonstrate_value_iteration(-1000)

# Case 3: Reward inside [-1, 1]
demonstrate_value_iteration(0.5)
