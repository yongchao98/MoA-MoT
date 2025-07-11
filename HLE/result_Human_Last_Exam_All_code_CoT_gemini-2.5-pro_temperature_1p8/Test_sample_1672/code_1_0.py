import numpy as np

def solve():
    """
    This function simulates the policy iteration process with an intrinsic reward
    to demonstrate the maximization of state entropy.
    """
    N_STATES = 10  # Number of states in a 1D grid
    
    # An initial policy that is heavily biased towards moving right
    policy = np.full(N_STATES, 0.95)

    max_entropy = np.log(N_STATES)
    print("This simulation demonstrates that the iterative policy update with an intrinsic reward")
    print(f"r_k(s) = -log(p_{{k-1}}(s)) leads to a policy that maximizes state entropy H(s).")
    print(f"The environment is a 1D grid with {N_STATES} states. The maximum possible entropy is log({N_STATES}) = {max_entropy:.4f}\n")
    print("{:<10} | {:<15} | {}".format("Iteration", "State Entropy", "Policy (Probability of going Right for each state)"))
    print("-" * 90)

    n_iterations = 15
    for k in range(n_iterations):
        # 1. For the current policy, calculate the state transition matrix T
        # T[i, j] = P(s_{t+1}=j | s_t=i)
        T = np.zeros((N_STATES, N_STATES))
        for s in range(N_STATES):
            p_right = policy[s]
            p_left = 1 - p_right
            
            # Action 'Right': move to s+1, or stay at boundary
            s_next_r = min(s + 1, N_STATES - 1)
            T[s, s_next_r] += p_right
            
            # Action 'Left': move to s-1, or stay at boundary
            s_next_l = max(s - 1, 0)
            T[s, s_next_l] += p_left

        # 2. Calculate the stationary state distribution p(s) for the current policy
        # p is the left eigenvector of T with eigenvalue 1. (p^T * T = p^T)
        # This is equivalent to the right eigenvector of T.T with eigenvalue 1.
        eigenvalues, eigenvectors = np.linalg.eig(T.T)
        idx = np.argmin(np.abs(eigenvalues - 1))
        p = np.real(eigenvectors[:, idx])
        p = p / np.sum(p) # Normalize to a probability distribution

        # 3. Calculate and print the entropy of the state distribution
        # Use a small epsilon to avoid log(0) for states with zero probability
        p_safe = p[p > 1e-12]
        entropy = -np.sum(p_safe * np.log(p_safe))
        policy_str = " ".join([f"{x:.2f}" for x in policy])
        print(f"{k:<10} | {entropy:<15.4f} | [{policy_str}]")

        # 4. Define the intrinsic reward r(s) = -log(p(s))
        rewards = -np.log(p + 1e-12)

        # 5. Update the policy for the next iteration.
        # The new policy pi_k will choose actions that lead to states with higher rewards.
        # We use a softmax function to determine the new probabilities.
        next_policy = np.zeros(N_STATES)
        for s in range(N_STATES):
            # One-step lookahead to get Q-values
            # Q-value for 'Right'
            s_next_r = min(s + 1, N_STATES - 1)
            q_right = rewards[s_next_r]
            
            # Q-value for 'Left'
            s_next_l = max(s - 1, 0)
            q_left = rewards[s_next_l]

            # Softmax to get new policy probability
            exp_q_right = np.exp(q_right)
            exp_q_left = np.exp(q_left)
            next_policy[s] = exp_q_right / (exp_q_right + exp_q_left)
            
        policy = next_policy

    print("-" * 90)
    print("As k increases, the state entropy converges to the maximum possible value.")
    print("The policy that achieves this is the limit of the iterative process, lim_{k->inf} pi^k.")


solve()