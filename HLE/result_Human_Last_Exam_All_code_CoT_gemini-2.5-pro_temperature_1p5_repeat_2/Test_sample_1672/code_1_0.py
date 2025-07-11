import numpy as np

def calculate_stationary_distribution(transition_matrix):
    """
    Calculates the stationary distribution of a Markov chain.
    The stationary distribution p is the left eigenvector of the transition matrix P
    with eigenvalue 1, i.e., pP = p. This is equivalent to finding the
    right eigenvector of P.T with eigenvalue 1.
    """
    try:
        eigenvalues, eigenvectors = np.linalg.eig(transition_matrix.T)
        # Find the eigenvector corresponding to the eigenvalue closest to 1
        stationary_vector_index = np.argmin(np.abs(eigenvalues - 1.0))
        stationary_vector = np.real(eigenvectors[:, stationary_vector_index])
        # Normalize the eigenvector to be a probability distribution
        stationary_distribution = stationary_vector / np.sum(stationary_vector)
        return stationary_distribution
    except np.linalg.LinAlgError:
        # Fallback for convergence issues: power iteration
        dist = np.ones(len(transition_matrix)) / len(transition_matrix)
        for _ in range(100):
            dist = dist.dot(transition_matrix)
        return dist


def calculate_entropy(probabilities):
    """Calculates the Shannon entropy of a probability distribution in bits."""
    # Filter out zero probabilities to avoid log(0)
    non_zero_probs = probabilities[probabilities > 1e-12]
    return -np.sum(non_zero_probs * np.log2(non_zero_probs))

def solve_max_entropy_policy_simulation():
    """
    Demonstrates that the policy converges to maximize state entropy through
    an iterative process using intrinsic motivation rewards.
    """
    # --- Environment Setup ---
    N_STATES = 7
    # Actions: 0 for 'left', 1 for 'right'
    N_ACTIONS = 2
    GAMMA = 0.99  # Discount factor for RL
    BETA = 10.0   # Inverse temperature for softmax policy improvement
    MAX_ITERATIONS = 15
    VALUE_ITERATION_TOLERANCE = 1e-6
    
    print("This simulation demonstrates the iterative process described in the problem.")
    print(f"The goal is to find a policy that maximizes the state entropy H(s).")
    print(f"The maximum possible entropy for {N_STATES} states is log2({N_STATES}) = {np.log2(N_STATES):.4f}\n")

    # --- Initial Policy (pi^0) ---
    # Start with a policy that heavily prefers going 'right'
    # policy_probs[s, a] = pi(a|s)
    policy_probs = np.zeros((N_STATES, N_ACTIONS))
    policy_probs[:, 1] = 1.0 # Go right with 100% probability

    # --- Main Iteration Loop ---
    for k in range(MAX_ITERATIONS):
        # The policy from the previous step is pi^k
        # We will analyze it and then compute the next policy pi^{k+1}
        
        # 1. Compute state transition matrix P_{pi^k}
        p_pi = np.zeros((N_STATES, N_STATES))
        for s in range(N_STATES):
            s_left = max(0, s - 1)
            s_right = min(N_STATES - 1, s + 1)
            p_pi[s, s_left] += policy_probs[s, 0]
            p_pi[s, s_right] += policy_probs[s, 1]
            
        # 2. Compute state distribution p_{pi^k}(s)
        state_distribution = calculate_stationary_distribution(p_pi)
        
        # 3. Calculate entropy H(p_{pi^k})
        entropy = calculate_entropy(state_distribution)
        print(f"Iteration k={k}:")
        print(f"  - State distribution for pi^{k}:   {np.round(state_distribution, 3)}")
        print(f"  - Entropy H(s) for pi^{k}:       {entropy:.4f}")

        # Now, we use this to find the policy for the next iteration, pi^{k+1}
        
        # 4. Define intrinsic reward r_{k+1}(s) = -log(p_{pi^k}(s))
        # Add a small epsilon to avoid log(0) when a state is unreachable
        rewards = -np.log(state_distribution + 1e-9)

        # 5. Find optimal Value function V*(s) for these rewards using Value Iteration
        V = np.zeros(N_STATES)
        while True:
            V_old = V.copy()
            Q = np.zeros((N_STATES, N_ACTIONS))
            for s in range(N_STATES):
                s_left = max(0, s - 1)
                s_right = min(N_STATES - 1, s + 1)
                # Reward is received based on the state you land in
                Q[s, 0] = rewards[s_left] + GAMMA * V_old[s_left]  # Q(s, 'left')
                Q[s, 1] = rewards[s_right] + GAMMA * V_old[s_right] # Q(s, 'right')
            V = np.max(Q, axis=1)
            if np.max(np.abs(V - V_old)) < VALUE_ITERATION_TOLERANCE:
                break
        
        # 6. Find new policy pi^{k+1} (Policy Improvement) via softmax
        final_Q = np.zeros((N_STATES, N_ACTIONS))
        for s in range(N_STATES):
            s_left = max(0, s - 1)
            s_right = min(N_STATES - 1, s + 1)
            final_Q[s, 0] = rewards[s_left] + GAMMA * V[s_left]
            final_Q[s, 1] = rewards[s_right] + GAMMA * V[s_right]
        exp_q = np.exp(BETA * (final_Q - np.max(final_Q, axis=1, keepdims=True)))
        policy_probs = exp_q / np.sum(exp_q, axis=1, keepdims=True)

    print("\nSimulation finished.")
    print("As k increases, the policy is updated to visit less frequent states.")
    print("This makes the state distribution p(s) more uniform, thus increasing the entropy H(s).")
    print("The process converges towards the maximum possible entropy.")
    print("\nConclusion: The policy that maximizes the entropy is the limit of this iterative process, lim_{k->inf} pi^k.")

solve_max_entropy_policy_simulation()
<<<A>>>