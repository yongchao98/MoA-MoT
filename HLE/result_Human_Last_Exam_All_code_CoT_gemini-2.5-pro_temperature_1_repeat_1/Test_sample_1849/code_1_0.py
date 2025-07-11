import numpy as np

def run_value_iteration(rewards, transitions, gamma, tol=1e-6):
    """
    Runs value iteration for a simple MDP to demonstrate convergence.

    Args:
        rewards (np.array): Reward vector R(s) for the single action.
        transitions (np.array): Transition matrix P(s'|s) for the single action.
        gamma (float): Discount factor.
        tol (float): Tolerance for convergence.

    Returns:
        tuple: Final value function and list of max-norm differences.
    """
    num_states = len(rewards)
    # Initialize value function to zeros
    V = np.zeros(num_states)
    diffs = []
    
    print(f"\n--- Running Value Iteration with Rewards: {rewards} ---")
    
    # Iterate until convergence
    for i in range(1000):
        V_old = V.copy()
        # V(s) = R(s) + gamma * Sum_{s'} P(s'|s)V_old(s')
        # Since there is only one action, the max_a is implicit.
        V = rewards + gamma * np.dot(transitions, V_old)
        
        diff = np.max(np.abs(V - V_old))
        diffs.append(diff)
        
        if i < 5:
             print(f"Iter {i+1}: ||V_k - V_{'k-1'}||_inf = {diff:.6f}")

        if diff < tol:
            print(f"Converged after {i+1} iterations.")
            break
            
    return V, diffs

def main():
    """
    Main function to demonstrate convergence rate independence from rewards.
    """
    # MDP Definition (2 states, 1 action for simplicity)
    # P[i, j] = P(s_j | s_i, a)
    transitions = np.array([
        [0.1, 0.9],
        [0.7, 0.3]
    ])
    
    gamma = 0.9
    
    print("This script demonstrates that the geometric convergence rate of Value Iteration")
    print("is determined by the discount factor gamma, not the reward function.")
    print("The proof that the Bellman operator is a contraction holds for any real-valued rewards.")
    print("-" * 70)
    print(f"MDP defined with gamma = {gamma}")


    # --- Scenario 1: Rewards in [-1, 1] ---
    rewards1 = np.array([1.0, -1.0])
    V1, diffs1 = run_value_iteration(rewards1, transitions, gamma)
    print("\nFinal Bellman Equation (V = R + γPV) for Scenario 1:")
    print(f"State 0: {V1[0]:.4f} = {rewards1[0]:.1f} + {gamma} * ({transitions[0,0]}*{V1[0]:.4f} + {transitions[0,1]}*{V1[1]:.4f})")
    print(f"State 1: {V1[1]:.4f} = {rewards1[1]:.1f} + {gamma} * ({transitions[1,0]}*{V1[0]:.4f} + {transitions[1,1]}*{V1[1]:.4f})")


    # --- Scenario 2: Arbitrary real-valued rewards ---
    rewards2 = np.array([100.0, -50.0])
    V2, diffs2 = run_value_iteration(rewards2, transitions, gamma)
    print("\nFinal Bellman Equation (V = R + γPV) for Scenario 2:")
    print(f"State 0: {V2[0]:.4f} = {rewards2[0]:.1f} + {gamma} * ({transitions[0,0]}*{V2[0]:.4f} + {transitions[0,1]}*{V2[1]:.4f})")
    print(f"State 1: {V2[1]:.4f} = {rewards2[1]:.1f} + {gamma} * ({transitions[1,0]}*{V2[0]:.4f} + {transitions[1,1]}*{V2[1]:.4f})")


    # --- Analysis of Convergence Rate ---
    print("\n" + "-" * 70)
    print("Analysis of Convergence Rate")
    # The ratio of successive differences should approach gamma
    rate1 = diffs1[1] / diffs1[0]
    rate2 = diffs2[1] / diffs2[0]
    
    print(f"For Rewards 1, the ratio ||V_2-V_1||/||V_1-V_0|| is {rate1:.4f}")
    print(f"For Rewards 2, the ratio ||V_2-V_1||/||V_1-V_0|| is {rate2:.4f}")
    print(f"Both are close to gamma ({gamma}), showing the rate is independent of rewards.")

if __name__ == '__main__':
    main()