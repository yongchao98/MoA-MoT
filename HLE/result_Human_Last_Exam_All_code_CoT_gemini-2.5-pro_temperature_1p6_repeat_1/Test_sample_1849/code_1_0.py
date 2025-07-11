import math

def demonstrate_value_iteration_convergence():
    """
    Demonstrates that value iteration convergence is independent of the reward scale.
    """
    # --- MDP Definition ---
    # States: S0, S1
    # In S0, only one action: 'go_to_s1', which deterministically lands in S1.
    # In S1, only one action: 'stay', which deterministically stays in S1.
    # This simplifies the Bellman update as the 'max_a' is over a single action.
    states = ["S0", "S1"]
    gamma = 0.9  # Discount factor
    theta = 0.001 # Convergence threshold

    def run_value_iteration(rewards, max_iterations=100):
        """
        Performs value iteration for the given rewards and prints each step.
        """
        print(f"--- Starting Value Iteration with Rewards: {rewards} ---")
        # Initialize value function
        V = {s: 0.0 for s in states}

        for i in range(max_iterations):
            delta = 0
            V_old = V.copy()
            
            print(f"Iteration {i+1}:")

            # Update for state S1
            # V(S1) = R(S1) + gamma * V_k(S1)
            old_v_s1 = V_old["S1"]
            new_v_s1 = rewards["S1"] + gamma * old_v_s1
            V["S1"] = new_v_s1
            print(f"  V(S1) update: V = {rewards['S1']} + {gamma} * {old_v_s1:.2f} = {V['S1']:.2f}")

            # Update for state S0
            # V(S0) = R(S0) + gamma * V_k(S1)
            old_v_s1_for_s0 = V_old["S1"] # Note: V(S0) depends on V(S1) from the previous iteration
            old_v_s0 = V_old["S0"]
            new_v_s0 = rewards["S0"] + gamma * old_v_s1_for_s0
            V["S0"] = new_v_s0
            print(f"  V(S0) update: V = {rewards['S0']} + {gamma} * {old_v_s1_for_s0:.2f} = {V['S0']:.2f}")

            # Check for convergence
            delta = max(abs(V["S0"] - old_v_s0), abs(V["S1"] - old_v_s1))
            print(f"  Max value change (delta): {delta:.4f}\n")
            
            if delta < theta:
                print(f"Convergence reached after {i+1} iterations.")
                break
        
        print(f"Final Value Function: S0={V['S0']:.2f}, S1={V['S1']:.2f}")
        print("--------------------------------------------------\n")

    # Case 1: Small rewards
    rewards_small = {"S0": 1, "S1": -1}
    run_value_iteration(rewards_small)

    # Case 2: Large rewards
    rewards_large = {"S0": 100, "S1": -100}
    run_value_iteration(rewards_large)

if __name__ == "__main__":
    demonstrate_value_iteration_convergence()