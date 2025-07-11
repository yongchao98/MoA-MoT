import numpy as np

def demonstrate_value_iteration_convergence():
    """
    Demonstrates that value iteration converges regardless of the reward range.
    """
    # 1. Define the Markov Decision Process (MDP)
    # A simple 2-state MDP
    num_states = 2
    num_actions = 2
    
    # Transition probabilities: P[state, action, next_state]
    # For simplicity, let's make transitions deterministic
    # Action 0 in state 0 -> state 0; Action 1 in state 0 -> state 1
    # Action 0 in state 1 -> state 1; Action 1 in state 1 -> state 0
    P = np.zeros((num_states, num_actions, num_states))
    P[0, 0, 0] = 1.0
    P[0, 1, 1] = 1.0
    P[1, 0, 1] = 1.0
    P[1, 1, 0] = 1.0

    # Reward function: R[state, action, next_state]
    # We use a wide range of rewards to demonstrate the principle.
    # Let's use rewards like -50 and +100.
    R = np.zeros((num_states, num_actions, num_states))
    R[0, 0, 0] = -50
    R[0, 1, 1] = 100
    R[1, 0, 1] = 100
    R[1, 1, 0] = -50
    
    # Discount factor (must be < 1 for guaranteed convergence)
    gamma = 0.9
    
    # Convergence threshold
    theta = 1e-6
    
    # 2. Initialize Value Function
    V = np.zeros(num_states)
    
    print("Starting Value Iteration...")
    print(f"MDP details: States={num_states}, Actions={num_actions}, Gamma={gamma}")
    print(f"Rewards are in the range [-50, 100].\n")

    iteration = 0
    while True:
        iteration += 1
        delta = 0
        V_old = V.copy()
        
        # 3. Perform Value Iteration Update for each state
        for s in range(num_states):
            # Calculate Q-values for each action
            q_values = np.zeros(num_actions)
            for a in range(num_actions):
                # Bellman equation: E[R + gamma * V(s')]
                q_values[a] = np.sum([P[s, a, s_prime] * (R[s, a, s_prime] + gamma * V_old[s_prime]) for s_prime in range(num_states)])
            
            # Update the value function with the max Q-value
            V[s] = np.max(q_values)
        
        # 4. Check for convergence
        delta = np.max(np.abs(V - V_old))
        print(f"Iteration {iteration}: Value Function = {np.round(V, 2)}, Change = {delta:.6f}")
        
        if delta < theta:
            break
            
    print("\nValue Iteration Converged!")
    print(f"Final Optimal Value Function: {V}")
    print("\nAs shown, convergence is achieved even with a wide range of rewards.")
    print("The convergence is guaranteed by the discount factor gamma < 1, not the reward values.")

demonstrate_value_iteration_convergence()