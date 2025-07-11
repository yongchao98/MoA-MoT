import numpy as np

def demonstrate_value_iteration_convergence():
    """
    Demonstrates that value iteration converges geometrically
    regardless of the scale of the reward function.
    """
    # MDP Parameters
    num_states = 2
    gamma = 0.9  # Discount factor
    # A simple transition model: from any state, 50/50 chance to move to the other state
    # This is a bit unrealistic but simple for demonstration. Let's make it simpler.
    # From state 0, go to state 1. From state 1, go to state 0.
    # T[s, s'] = P(s' | s, only_action)
    T = np.array([[0.0, 1.0], [1.0, 0.0]])
    
    # Rewards R(s). We use values outside [-1, 1] to prove the point.
    R = np.array([100.0, -50.0])

    print(f"Demonstrating convergence with large rewards R = {R.tolist()} and gamma = {gamma}\n")
    print("Value Iteration Steps:")
    print("----------------------")
    print(f"{'Iteration':<10} | {'V(s0)':<15} | {'V(s1)':<15} | {'Max Difference':<20}")
    print("-" * 65)

    # Initialize value function to zeros
    V = np.zeros(num_states)
    
    for i in range(20):
        V_old = V.copy()
        
        # In this simple model, there's only one action from each state, so no 'max' is needed.
        # V_new(s) = R(s) + gamma * sum(P(s'|s) * V_old(s'))
        # This can be written in vector form
        V = R + gamma * T @ V_old

        # Calculate the max norm of the difference
        max_diff = np.max(np.abs(V - V_old))
        
        # Print the values and the difference
        print(f"{i+1:<10} | {V[0]:<15.4f} | {V[1]:<15.4f} | {max_diff:<20.4f}")

        # The difference should approach 0
        if max_diff < 1e-4:
            print("\nConvergence reached.")
            break
            
    # Calculate the theoretical optimal value function V*
    # V* = (I - gamma*T)^-1 * R
    I = np.identity(num_states)
    V_star = np.linalg.inv(I - gamma * T) @ R
    print(f"\nTheoretical optimal value function V*: {V_star.tolist()}")
    print(f"Final iterated value function V:      {V.tolist()}")

demonstrate_value_iteration_convergence()