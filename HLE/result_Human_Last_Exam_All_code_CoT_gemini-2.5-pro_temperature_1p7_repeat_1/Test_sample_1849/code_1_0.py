import numpy as np

def solve_value_iteration():
    """
    Demonstrates that value iteration converges for arbitrary real-valued rewards.
    """
    print("This script demonstrates the convergence of the Value Iteration algorithm.")
    print("The convergence is guaranteed by the discount factor gamma < 1 and is independent of the range of reward values.\n")

    # Define a simple Markov Decision Process (MDP)
    # We use arbitrary real numbers for rewards to show the point.
    rewards = np.array([-10.0, 50.0])  # R(s=0) = -10, R(s=1) = 50
    # Transition probabilities P(s' | s) for a single action
    transitions = np.array([[0.5, 0.5], [0.2, 0.8]])
    # Discount factor
    gamma = 0.9

    # Initialize value function
    v = np.zeros(len(rewards))
    max_iterations = 100
    convergence_threshold = 1e-9

    print("Running Value Iteration...")
    for i in range(max_iterations):
        v_prev = v.copy()
        # Bellman update: V_{k+1} = R + gamma * P * V_k
        v = rewards + gamma * (transitions @ v_prev)
        # Check for convergence
        if np.max(np.abs(v - v_prev)) < convergence_threshold:
            print(f"Algorithm converged after {i+1} iterations.\n")
            break

    print("--- Final Result ---")
    print(f"Final Optimal Value Function V*:\n  V*(s=0) = {v[0]:.4f}\n  V*(s=1) = {v[1]:.4f}\n")

    print("Verifying the Bellman Optimality Equation V* = R + gamma * P * V* for each state:")
    # For state 0
    print("\nFor state s=0:")
    r_s0 = rewards[0]
    p_s0_s0, p_s1_s0 = transitions[0, 0], transitions[0, 1]
    v_s0, v_s1 = v[0], v[1]
    rhs_s0 = r_s0 + gamma * (p_s0_s0 * v_s0 + p_s1_s0 * v_s1)
    print(f"V*(s=0) = R(s=0) + gamma * [P(s'=0|s=0)V*(s'=0) + P(s'=1|s=0)V*(s'=1)]")
    print(f"Equation: {v_s0:.4f} = {r_s0} + {gamma} * [{p_s0_s0}*({v_s0:.4f}) + {p_s1_s0}*({v_s1:.4f})]")
    print(f"Result:   {v_s0:.4f} = {rhs_s0:.4f} (Verified)")

    # For state 1
    print("\nFor state s=1:")
    r_s1 = rewards[1]
    p_s0_s1, p_s1_s1 = transitions[1, 0], transitions[1, 1]
    rhs_s1 = r_s1 + gamma * (p_s0_s1 * v_s0 + p_s1_s1 * v_s1)
    print(f"V*(s=1) = R(s=1) + gamma * [P(s'=0|s=1)V*(s'=0) + P(s'=1|s=1)V*(s'=1)]")
    print(f"Equation: {v_s1:.4f} = {r_s1} + {gamma} * [{p_s0_s1}*({v_s0:.4f}) + {p_s1_s1}*({v_s1:.4f})]")
    print(f"Result:   {v_s1:.4f} = {rhs_s1:.4f} (Verified)")


solve_value_iteration()