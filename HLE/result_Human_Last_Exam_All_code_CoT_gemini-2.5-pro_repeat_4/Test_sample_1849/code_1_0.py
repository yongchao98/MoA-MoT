import numpy as np

def demonstrate_value_iteration_convergence():
    """
    Demonstrates that value iteration converges even with large,
    unbounded-seeming rewards in a finite state-space MDP.
    """
    # 1. Define the MDP (S={s0, s1}, A={a0, a1})
    num_states = 2
    num_actions = 2
    gamma = 0.9  # Discount factor

    # Rewards R(s, a)
    # We use large rewards to show convergence is independent of their scale.
    R = np.array([
        [1.0, 1000.0],     # R(s0, a0), R(s0, a1)
        [2.0, -5000.0]     # R(s1, a0), R(s1, a1)
    ])

    # Transition probabilities P(s'|s, a) are deterministic for simplicity
    # P is a |S|x|A|x|S| matrix
    P = np.zeros((num_states, num_actions, num_states))
    P[0, 0, 0] = 1.0  # s0, a0 -> s0
    P[0, 1, 1] = 1.0  # s0, a1 -> s1
    P[1, 0, 1] = 1.0  # s1, a0 -> s1
    P[1, 1, 0] = 1.0  # s1, a1 -> s0

    # 2. Value Iteration Algorithm
    V = np.zeros(num_states) # Initialize V(s) = 0 for all s
    tolerance = 1e-6
    max_iterations = 1000

    print("Running Value Iteration...")
    for i in range(max_iterations):
        V_old = V.copy()
        Q = np.zeros((num_states, num_actions))
        for s in range(num_states):
            for a in range(num_actions):
                # E[V(s')] = Σ_{s'} P(s'|s,a)V(s')
                expected_future_value = P[s, a, :].dot(V_old)
                Q[s, a] = R[s, a] + gamma * expected_future_value

        V = np.max(Q, axis=1)

        # Check for convergence
        if np.max(np.abs(V - V_old)) < tolerance:
            print(f"Converged after {i+1} iterations.\n")
            break

    # 3. Output the final result
    V_star = V
    print("The algorithm converges to the optimal value function V*:")
    for s in range(num_states):
        print(f"V*(s{s}) = {V_star[s]:.4f}")

    print("\nVerifying the Bellman Optimality Equation with the final values:")
    for s in range(num_states):
        print(f"\nFor state s{s}:")
        # For s0, the optimal action is a1. For s1, it is a0.
        # R(s0,a0) + g*V*(s0), R(s0,a1) + g*V*(s1)
        val_a0 = R[s,0] + gamma * P[s, 0, :].dot(V_star)
        # R(s1,a0) + g*V*(s1), R(s1,a1) + g*V*(s0)
        val_a1 = R[s,1] + gamma * P[s, 1, :].dot(V_star)
        
        print(f"V*(s{s}) = {V_star[s]:.4f}")
        print(f"  max[ R(s{s},a0) + γ*V*(s'), R(s{s},a1) + γ*V*(s') ]")
        print(f"  max[ {R[s,0]} + {gamma}*V_next_s_a0, {R[s,1]} + {gamma}*V_next_s_a1 ]")
        if s == 0:
            print(f"  max[ {R[s,0]} + {gamma}*{V_star[0]:.4f}, {R[s,1]} + {gamma}*{V_star[1]:.4f} ]")
        else: # s == 1
            print(f"  max[ {R[s,0]} + {gamma}*{V_star[1]:.4f}, {R[s,1]} + {gamma}*{V_star[0]:.4f} ]")

        print(f"  max[ {val_a0:.4f}, {val_a1:.4f} ] = {np.max([val_a0, val_a1]):.4f}")


demonstrate_value_iteration_convergence()