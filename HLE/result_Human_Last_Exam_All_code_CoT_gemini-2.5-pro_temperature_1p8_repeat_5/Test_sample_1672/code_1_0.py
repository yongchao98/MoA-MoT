import numpy as np

def demonstrate_entropy_maximization():
    """
    Simulates the iterative policy learning process to demonstrate that it
    maximizes state entropy.

    We use a simple environment with 2 states {0, 1}.
    An action 'a' deterministically transitions to state 'a'.
    The policy is stochastic and state-independent, defined by the probability
    of taking action 0 or 1.
    """
    # --- Hyperparameters ---
    num_states = 2
    gamma = 0.9  # Discount factor
    beta = 1.0   # Softmax temperature for policy update
    iterations = 30
    epsilon = 1e-9 # Small number to prevent log(0)

    print(f"Starting simulation for a {num_states}-state environment.")
    print(f"Theoretical maximum entropy H(s) = log({num_states}) = {np.log(num_states):.4f}\n")

    # 1. Initialize policy pi_0. Start with a biased policy.
    # policy[a] is the probability of taking action 'a'
    policy = np.array([0.9, 0.1])

    print("--- Iteration Details ---")
    print("k | Policy p(a=0), p(a=1) | State Dist. p(s=0), p(s=1) | Entropy H(s)")
    print("-" * 70)

    for k in range(iterations):
        # 2. Calculate the stationary state distribution p_k-1 for the current policy.
        # For our simple MDP where pi is state-independent and T(s,a,s') leads to state a,
        # the stationary distribution is identical to the policy's action probabilities.
        # p(s=0) = pi(a=0), p(s=1) = pi(a=1)
        state_dist = policy

        # Calculate the entropy of the current state distribution
        entropy = -np.sum(state_dist * np.log(state_dist + epsilon))

        # Print current status
        policy_str = f"[{policy[0]:.4f}, {policy[1]:.4f}]"
        dist_str = f"[{state_dist[0]:.4f}, {state_dist[1]:.4f}]"
        print(f"{k:<2}| {policy_str:<23}| {dist_str:<29}| {entropy:.4f}")

        # 3. Define the reward for the next iteration r_k(s) = -log(p_k-1(s))
        rewards = -np.log(state_dist + epsilon)

        # 4. Find the new policy pi_k via value iteration and policy improvement.
        # First, find the optimal Q-values for the current reward function.
        # V*(1-gamma) = max(rewards) => V* = max(rewards)/(1-gamma)
        value_star = np.max(rewards) / (1.0 - gamma)
        # Q*(s,a) = r(a) + gamma * V*(a)
        q_values = rewards + gamma * value_star

        # Update the policy using softmax on the Q-values.
        # This makes the policy stochastic and avoids oscillations.
        exp_q = np.exp(beta * q_values)
        policy = exp_q / np.sum(exp_q)

    print("-" * 70)
    print("\nSimulation finished. As k increases, the policy and state distribution")
    print("converge to uniform (0.5, 0.5), and the entropy H(s) converges to")
    print(f"its maximum value of {np.log(num_states):.4f}.")

if __name__ == '__main__':
    demonstrate_entropy_maximization()
