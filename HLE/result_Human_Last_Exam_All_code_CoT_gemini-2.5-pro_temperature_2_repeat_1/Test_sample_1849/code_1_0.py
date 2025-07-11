import numpy as np

def value_iteration_demo(R, gamma=0.9, T=None, max_iterations=100, tolerance=1e-4):
    """
    Demonstrates value iteration on a simple 2-state, 2-action MDP.

    Args:
        R (np.array): Reward matrix of shape (actions, states).
        gamma (float): Discount factor.
        T (np.array, optional): Transition matrix of shape (actions, start_states, end_states).
        max_iterations (int): Maximum number of iterations.
        tolerance (float): Convergence threshold.
    """
    if T is None:
        # T[action, start_state, end_state]
        # Action 0: tend to stay, Action 1: tend to move
        T = np.array([
            [[0.9, 0.1], [0.1, 0.9]],  # Action 0
            [[0.1, 0.9], [0.9, 0.1]]   # Action 1
        ])

    num_states = T.shape[1]
    V = np.zeros(num_states)
    
    print(f"--- Running Value Iteration ---")
    print(f"Rewards R[action, state]:\n{R}")
    print("---------------------------------------------")

    for i in range(max_iterations):
        V_old = V.copy()
        
        # Calculate Q-values: Q(s,a) = R(s,a) + gamma * Sum[P(s'|s,a) * V_old(s')]
        # Using numpy broadcasting and einsum for efficiency
        Q = R + gamma * np.einsum('ijk,k->ij', T, V_old)
        
        # Update value function: V(s) = max_a Q(s,a)
        V = np.max(Q, axis=0)
        
        # Check for convergence
        delta = np.max(np.abs(V - V_old))
        print(f"Iteration {i+1:3}: V = {np.array2string(V, precision=4, sign=' ')}, Delta = {delta:.6f}")
        
        # Stop if converged or if values explode
        if np.isinf(V).any():
            print("\nValues have become infinite. Divergence.")
            break
        if delta < tolerance:
            print(f"\nConverged to a stable value function after {i+1} iterations.")
            break
    else:
        # This part runs if the loop finishes without break
        print(f"\nStopped after {max_iterations} iterations without converging to the tolerance.")
        
    print(f"Final Value Function: V = {np.array2string(V, precision=4, sign=' ')}")
    print("-" * 45)


# Case 1: Bounded rewards (within the range [-1, 1])
print("DEMONSTRATION 1: Bounded Rewards from [-1, 1]\n")
R_bounded = np.array([
    [-1.0,  1.0],  # Rewards for action 0 in states 0 and 1
    [ 0.5, -0.5]   # Rewards for action 1 in states 0 and 1
])
value_iteration_demo(R_bounded)

# Case 2: Unbounded rewards
print("\n\nDEMONSTRATION 2: Unbounded Rewards (from R)\n")
R_unbounded = np.array([
    [-1.0, float('inf')], # Infinite reward for (action=0, state=1)
    [ 0.5, -0.5]
])
# We run for fewer iterations just to see the values explode to infinity
value_iteration_demo(R_unbounded, max_iterations=10)
