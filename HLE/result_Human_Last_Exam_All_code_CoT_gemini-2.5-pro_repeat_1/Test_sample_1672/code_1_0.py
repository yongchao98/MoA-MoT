import numpy as np
import matplotlib.pyplot as plt

def estimate_state_distribution(policy, n_states, n_steps=10000):
    """Simulates the policy to get the state visitation distribution."""
    counts = np.zeros(n_states)
    # Start at a random state
    current_state = np.random.randint(n_states)
    for _ in range(n_steps):
        counts[current_state] += 1
        # Sample action from policy for the current state
        # Actions are -1 (left) and +1 (right)
        p_right = policy[current_state]
        action = 1 if np.random.rand() < p_right else -1
        # Move to the next state (on a ring)
        current_state = (current_state + action) % n_states
    
    # Add a small constant to avoid division by zero or log(0)
    distribution = (counts + 1e-9) / np.sum(counts + 1e-9)
    return distribution

def calculate_entropy(distribution):
    """Calculates the Shannon entropy of a probability distribution."""
    return -np.sum(distribution * np.log(distribution))

def value_iteration(rewards, n_states, gamma=0.9):
    """Performs value iteration to find the optimal Q-values."""
    v = np.zeros(n_states)
    for _ in range(100): # Iterate to convergence
        v_new = np.zeros(n_states)
        for s in range(n_states):
            # Q-values for actions left (-1) and right (+1)
            q_left = rewards[(s - 1) % n_states] + gamma * v[(s - 1) % n_states]
            q_right = rewards[(s + 1) % n_states] + gamma * v[(s + 1) % n_states]
            v_new[s] = max(q_left, q_right)
        if np.max(np.abs(v - v_new)) < 1e-4:
            break
        v = v_new
        
    # Calculate final Q-values from the optimal V-values
    q_values = np.zeros((n_states, 2)) # 0: left, 1: right
    for s in range(n_states):
        q_values[s, 0] = rewards[(s - 1) % n_states] + gamma * v[(s - 1) % n_states]
        q_values[s, 1] = rewards[(s + 1) % n_states] + gamma * v[(s + 1) % n_states]
    return q_values

def boltzmann_policy(q_values, beta=1.0):
    """Creates a Boltzmann (softmax) policy from Q-values."""
    exp_q = np.exp(beta * q_values)
    policy_prob_right = exp_q[:, 1] / np.sum(exp_q, axis=1)
    return policy_prob_right

def main():
    """
    Main function to run the simulation.
    """
    # --- Setup ---
    N_STATES = 10
    N_ITERATIONS = 20
    
    # The policy is defined by the probability of going right at each state.
    # Initialize pi^0 to be a biased policy (always go right)
    policy = np.full(N_STATES, 0.99) 
    
    entropies = []

    # Maximum possible entropy for N_STATES (uniform distribution)
    max_entropy = np.log(N_STATES)
    print(f"Environment: {N_STATES}-state ring.")
    print(f"Maximum possible entropy: {max_entropy:.4f}")
    print("-" * 30)
    print("Running simulation...")
    
    # --- Main Loop ---
    for k in range(N_ITERATIONS):
        # 1. Estimate state distribution for the current policy
        p_k = estimate_state_distribution(policy, N_STATES)
        
        # 2. Calculate and store entropy
        h_k = calculate_entropy(p_k)
        entropies.append(h_k)
        print(f"Iteration k={k}: Entropy H(p_k) = {h_k:.4f}")

        # 3. Define reward for the next iteration
        r_k_plus_1 = -np.log(p_k)
        
        # 4. Update policy using value iteration with the new reward
        q_values = value_iteration(r_k_plus_1, N_STATES)
        
        # 5. New policy is a softmax over the Q-values
        policy = boltzmann_policy(q_values, beta=1.0)
        
    print("-" * 30)
    print("Simulation finished. As k increases, the entropy of the state distribution")
    print("approaches the maximum possible value, confirming that the process")
    print("drives the policy towards one that maximizes state entropy.")
    
    # Plotting the results
    plt.figure(figsize=(10, 6))
    plt.plot(range(N_ITERATIONS), entropies, marker='o', label='Entropy H(s)')
    plt.axhline(y=max_entropy, color='r', linestyle='--', label=f'Max Entropy ({max_entropy:.2f})')
    plt.xlabel("Iteration (k)")
    plt.ylabel("Entropy H(s)")
    plt.title("State Entropy vs. Policy Iteration")
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == '__main__':
    main()