import numpy as np

def demonstrate_value_iteration_convergence():
    """
    Demonstrates that value iteration converges geometrically regardless
    of the range of the reward function.
    """
    # 1. Define the MDP components
    # Two states: S0, S1
    # Two actions: A0, A1
    num_states = 2
    num_actions = 2

    # Transition probabilities: P[s, a, s'] = P(s' | s, a)
    # A simplified model where actions have deterministic outcomes for clarity
    P = np.zeros((num_states, num_actions, num_states))
    P[0, 0, 0] = 1.0  # In S0, taking A0 stays in S0
    P[0, 1, 1] = 1.0  # In S0, taking A1 goes to S1
    P[1, 0, 1] = 1.0  # In S1, taking A0 stays in S1
    P[1, 1, 0] = 1.0  # In S1, taking A1 goes to S0

    # Reward function: R[s, a]
    # We use large positive and negative rewards to show it doesn't affect convergence.
    R = np.array([
        [1000, -500],   # Rewards in S0 for A0, A1
        [-2000, 800]    # Rewards in S1 for A0, A1
    ])
    
    # Discount factor (must be < 1 for convergence)
    gamma = 0.9

    # 2. Value Iteration Algorithm
    # Initialize value function to zeros
    V = np.zeros(num_states)
    
    print(f"Starting value iteration with gamma = {gamma}")
    print(f"Reward function R:\n{R}\n")
    print("Iteration | Max Value-Function Change")
    print("--------------------------------------")

    max_iterations = 20
    for i in range(max_iterations):
        V_old = V.copy()
        Q = np.zeros((num_states, num_actions))
        
        # Calculate the Q-values for all state-action pairs
        for s in range(num_states):
            for a in range(num_actions):
                # Expected future value: Σ_{s'} P(s' | s, a) * V_old(s')
                expected_future_value = P[s, a, :].dot(V_old)
                Q[s, a] = R[s, a] + gamma * expected_future_value
        
        # Update the value function by taking the max over actions
        V = np.max(Q, axis=1)
        
        # Calculate the maximum change (using the infinity norm)
        max_change = np.max(np.abs(V - V_old))
        
        print(f"{i+1:9d} | {max_change:.6f}")

        # Check for convergence
        if max_change < 1e-6:
            print("\nConvergence reached.")
            break
    
    print("\nFinal Optimal Value Function:")
    print(f"V*(S0) = {V[0]:.2f}")
    print(f"V*(S1) = {V[1]:.2f}")
    print("\nNotice how the max change decreases geometrically, proving convergence")
    print("is independent of the reward values' range.")


demonstrate_value_iteration_convergence()