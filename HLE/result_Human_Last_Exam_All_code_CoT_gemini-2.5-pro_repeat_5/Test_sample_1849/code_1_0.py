import numpy as np

def demonstrate_value_iteration_convergence(rewards, gamma, tolerance=1e-5):
    """
    Demonstrates the convergence of value iteration for a given reward function.

    Args:
        rewards (list): A list of rewards for each state.
        gamma (float): The discount factor.
        tolerance (float): The convergence threshold.
    """
    print(f"\n--- Running experiment with rewards: {rewards} ---")

    # Simple MDP setup:
    # 2 states (S0, S1)
    # 1 action (implicit, since there's only one choice)
    # Transition probabilities P(s' | s)
    P = np.array([
        [0.7, 0.3],  # Transitions from S0
        [0.2, 0.8]   # Transitions from S1
    ])
    R = np.array(rewards)

    # Initialize value function
    V = np.zeros(len(R))
    print(f"Initial V: {V}")

    iteration = 0
    while True:
        iteration += 1
        V_old = V.copy()

        # Bellman update for each state
        # V_new[s] = R[s] + gamma * sum(P[s, s'] * V_old[s'])
        expected_future_value = P @ V_old
        V = R + gamma * expected_future_value

        # Print the update equation for transparency
        print(f"\nIteration {iteration}:")
        for s in range(len(R)):
            # This loop builds the string showing the calculation for each state
            future_val_calc_terms = []
            for s_prime in range(len(R)):
                term = f"{P[s, s_prime]:.1f} * {V_old[s_prime]:.2f}"
                future_val_calc_terms.append(term)
            future_val_calc = " + ".join(future_val_calc_terms)
            
            print(f"  V_new(S{s}) = {R[s]} + {gamma} * ({future_val_calc}) = {V[s]:.4f}")

        # Check for convergence by looking at the max-norm difference
        max_diff = np.max(np.abs(V - V_old))
        print(f"-> Max difference ||V_new - V_old||_inf = {max_diff:.6f}")

        if max_diff < tolerance:
            print(f"\nConverged after {iteration} iterations.")
            print(f"Final Value Function: {V}")
            break
        
        if iteration > 50: # safety break to prevent infinite loops
            print("Did not converge in 50 iterations.")
            break

# Run the demonstration with two different reward scales.
# The discount factor gamma is 0.9.
demonstrate_value_iteration_convergence(rewards=[1, -1], gamma=0.9)
demonstrate_value_iteration_convergence(rewards=[1000, -1000], gamma=0.9)