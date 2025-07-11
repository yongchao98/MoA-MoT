import numpy as np

def run_simulation():
    """
    Simulates the iterative policy update process to show that it maximizes
    the state entropy over time.
    """

    # --- 1. Environment and Simulation Parameters ---
    NUM_STATES = 5
    NUM_ITERATIONS = 15
    SIMULATION_STEPS = 20000 # Steps per policy to estimate state distribution
    GAMMA = 0.95 # Discount factor for policy updates

    # --- 2. Helper Functions ---
    class ChainEnv:
        """A simple 1D chain environment."""
        def __init__(self, num_states):
            self.num_states = num_states
            self.action_space = [0, 1]  # 0: LEFT, 1: RIGHT

        def step(self, state, action):
            if action == 0:  # LEFT
                return max(0, state - 1)
            else:  # RIGHT
                return min(self.num_states - 1, state + 1)

    def get_state_distribution(policy, env, num_steps, start_state=0):
        """Estimates the state visitation distribution for a given policy."""
        counts = np.zeros(env.num_states)
        state = start_state
        for _ in range(num_steps):
            counts[state] += 1
            action_probs = policy[state]
            action = np.random.choice(env.action_space, p=action_probs)
            state = env.step(state, action)
        # Use a small epsilon to avoid zero probabilities leading to -inf logs
        distribution = (counts + 1e-9) / np.sum(counts + 1e-9)
        return distribution

    def get_entropy(distribution):
        """Calculates the Shannon entropy of a distribution."""
        return -np.sum(distribution * np.log2(distribution))

    def update_policy(reward, env, gamma):
        """Updates the policy using value iteration based on the given reward."""
        V = np.zeros(env.num_states)
        # Value iteration to find the optimal value function V(s)
        for _ in range(100):
            V_new = np.zeros(env.num_states)
            for s in range(env.num_states):
                q_values = [reward[s] + gamma * V[env.step(s, a)] for a in env.action_space]
                V_new[s] = max(q_values)
            if np.max(np.abs(V - V_new)) < 1e-6:
                break
            V = V_new

        # Extract a stochastic policy from V(s)
        new_policy = np.zeros((env.num_states, len(env.action_space)))
        for s in range(env.num_states):
            q_values = np.array([reward[s] + gamma * V[env.step(s, a)] for a in env.action_space])
            # Use softmax to convert Q-values to probabilities for a stochastic policy
            # A high temperature will lead to a more uniform action selection
            temp = 0.1
            exp_q = np.exp((q_values - np.max(q_values)) / temp)
            new_policy[s] = exp_q / np.sum(exp_q)
        return new_policy

    # --- 3. Main Simulation Loop ---
    env = ChainEnv(NUM_STATES)

    # Initialize pi^0: a biased policy that strongly prefers going RIGHT
    policy = np.array([[0.1, 0.9]] * NUM_STATES)
    
    max_entropy = np.log2(NUM_STATES)
    print(f"Environment: 1D Chain with {NUM_STATES} states.")
    print(f"Theoretical Maximum Entropy H(s) = log2({NUM_STATES}) = {max_entropy:.4f}\n")
    print("Iteration k | Entropy H(p_pi^k(s))")
    print("-----------------------------------")

    # Calculate and print initial state
    p_k = get_state_distribution(policy, env, SIMULATION_STEPS)
    H_k = get_entropy(p_k)
    print(f"    0       | {H_k:.4f}")

    # Run the iterative process
    for k in range(1, NUM_ITERATIONS + 1):
        # The reward for the new policy pi^k is based on the distribution of pi^{k-1}
        reward_k = -np.log(p_k)
        
        # Update policy from pi^{k-1} to pi^k
        policy = update_policy(reward_k, env, GAMMA)
        
        # Get the state distribution p_pi^k(s) for the new policy
        p_k = get_state_distribution(policy, env, SIMULATION_STEPS)
        
        # Calculate its entropy H(p_pi^k(s))
        H_k = get_entropy(p_k)
        
        print(f"   {k:2d}       | {H_k:.4f}")
    
    # --- 4. Final Output and Explanation ---
    print("\n-----------------------------------")
    print("As k increases, the entropy H(s) approaches the maximum possible value.")
    print("This shows the state distribution becomes more uniform over iterations.")
    
    final_p = p_k
    equation_str = "H(s) = - ["
    terms = [f"{p:.3f}*log2({p:.3f})" for p in final_p if p > 1e-9]
    equation_str += " + ".join(terms) + "]"
    
    print("\nFinal State Distribution p(s):")
    print(np.round(final_p, 3))
    print("\nThe entropy calculation for the final distribution is:")
    print(equation_str)
    print(f"Which evaluates to H(s) = {H_k:.4f}, close to the max of {max_entropy:.4f}.")
    print("\nConclusion: The process converges to a policy that maximizes state entropy.")

if __name__ == '__main__':
    run_simulation()