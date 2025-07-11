import numpy as np

def get_state_distribution(policy, num_states, num_steps=50000):
    """Simulates the policy to get the state distribution."""
    counts = np.zeros(num_states)
    state = np.random.randint(num_states) # Start at a random state
    for _ in range(num_steps):
        counts[state] += 1
        if np.random.rand() < policy[state]['right']:
            state = min(state + 1, num_states - 1)
        else:
            state = max(state - 1, 0)
    # Add a small epsilon to avoid division by zero or log(0)
    p = (counts + 1e-9) / np.sum(counts + 1e-9)
    return p

def calculate_entropy(p):
    """Calculates the entropy of a distribution using log base 2."""
    return -np.sum(p * np.log2(p))

def main():
    """
    Simulates the policy iteration process to show entropy maximization.
    """
    # --- Parameters ---
    NUM_STATES = 10
    NUM_ITERATIONS = 15

    # --- Initialization (k=0) ---
    # Start with a policy pi^0 that strongly prefers going right.
    # policy[s] = {'right': prob_right, 'left': prob_left}
    policy_k = [{'right': 0.95, 'left': 0.05} for _ in range(NUM_STATES)]
    # At the rightmost state, it must go left to avoid getting stuck.
    policy_k[NUM_STATES - 1] = {'right': 0.0, 'left': 1.0}
    
    print(f"Goal: Maximize state entropy H(s).")
    print(f"A uniform distribution over {NUM_STATES} states has max entropy H(s) = log2({NUM_STATES}) = {np.log2(NUM_STATES):.4f}")
    print("-" * 60)
    print("This simulation shows how H(s) increases with each policy iteration k.")
    print("-" * 60)
    
    # Calculate initial entropy for policy pi^0
    p_k = get_state_distribution(policy_k, NUM_STATES)
    entropy_k = calculate_entropy(p_k)
    print(f"Iteration k=0: H(s) = {entropy_k:.4f}")
    
    # --- Iteration Loop ---
    # `p_k` is p^{k-1} from the previous iteration
    for k in range(1, NUM_ITERATIONS + 1):
        # 1. Current state distribution is p_k (from pi^{k-1})
        p_prev = p_k

        # 2. Calculate intrinsic rewards r_k(s) = -log(p^{k-1}(s))
        # High reward for states that were rarely visited
        rewards = -np.log(p_prev)
        
        # 3. Update policy to pi^k.
        # The new policy will favor actions leading to states with higher rewards.
        # We model this using a softmax over the rewards of the potential next states.
        policy_next = [{} for _ in range(NUM_STATES)]
        for s in range(NUM_STATES):
            reward_left = rewards[max(0, s-1)]
            reward_right = rewards[min(NUM_STATES-1, s+1)]
            
            # Use softmax to convert rewards to action probabilities
            # This makes the policy stochastic, preferring high-reward next states
            exp_r = np.exp(reward_right)
            exp_l = np.exp(reward_left)
            prob_right = exp_r / (exp_r + exp_l)
            
            # Handle boundary conditions
            if s == 0:
                prob_right = 1.0 # Must go right from the leftmost state
            elif s == NUM_STATES - 1:
                prob_right = 0.0 # Must go left from the rightmost state

            policy_next[s] = {'right': prob_right, 'left': 1.0 - prob_right}
        
        policy_k = policy_next

        # 4. Calculate new state distribution p^k from pi^k and its entropy H(s)
        p_k = get_state_distribution(policy_k, NUM_STATES)
        entropy_k = calculate_entropy(p_k)

        print(f"Iteration k={k}: H(s) = {entropy_k:.4f}")

if __name__ == '__main__':
    main()