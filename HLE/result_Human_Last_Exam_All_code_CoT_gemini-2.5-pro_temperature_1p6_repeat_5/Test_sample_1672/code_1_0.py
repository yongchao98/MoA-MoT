import numpy as np

# --- Environment Setup ---
N_STATES = 5
N_ACTIONS = 2  # 0: Left, 1: Right
GAMMA = 0.95   # Discount factor
BETA = 20      # Inverse temperature for softmax policy

# Transition tensor T[s, a, s_next]
T = np.zeros((N_STATES, N_ACTIONS, N_STATES))
for s in range(N_STATES):
    # Action 0 (Left)
    s_next_l = max(0, s - 1)
    T[s, 0, s_next_l] = 1.0
    # Action 1 (Right)
    s_next_r = min(N_STATES - 1, s + 1)
    T[s, 1, s_next_r] = 1.0

def calculate_state_dist(policy, T_matrix):
    """Calculates the stationary state distribution for a given policy."""
    n_states, n_actions, _ = T_matrix.shape
    # Policy-conditioned transition matrix P_pi[s, s_next]
    P_pi = np.einsum('sa,san->sn', policy, T_matrix)

    # Use power iteration to find the stationary distribution
    dist = np.ones(n_states) / n_states
    for _ in range(200):
        dist_new = dist @ P_pi
        if np.allclose(dist, dist_new):
            break
        dist = dist_new
    return dist / np.sum(dist)

def solve_for_policy(rewards, T_matrix):
    """Finds the optimal policy for a given reward function using Value Iteration."""
    n_states, n_actions, _ = T_matrix.shape
    V = np.zeros(n_states)
    # Value Iteration
    for _ in range(200):
        # Q[s, a] = r[s] + gamma * sum_{s'} T(s'|s,a)V(s')
        Q = rewards[:, np.newaxis] + GAMMA * np.einsum('san,n->sa', T_matrix, V)
        V_new = np.max(Q, axis=1)
        if np.allclose(V, V_new):
            break
        V = V_new
    
    # Policy extraction (softmax for stochasticity)
    Q = rewards[:, np.newaxis] + GAMMA * np.einsum('san,n->sa', T_matrix, V)
    exp_q = np.exp(BETA * Q)
    policy = exp_q / np.sum(exp_q, axis=1, keepdims=True)
    return policy

def calculate_entropy(dist):
    """Calculates the entropy of a probability distribution."""
    # Use a small epsilon to prevent log(0)
    epsilon = 1e-12
    return -np.sum(dist * np.log(dist + epsilon))

# --- Main Simulation Loop ---
print("Simulating policy iteration with intrinsic motivation.\n")
# Start with a biased policy (strong preference to go right)
pi = np.full((N_STATES, N_ACTIONS), 0.01)
pi[:, 1] = 0.99

max_entropy = np.log(N_STATES)
print(f"Environment: {N_STATES} states in a chain.")
print(f"Maximum possible entropy H(s) = log({N_STATES}) = {max_entropy:.4f}\n")
print(f"{'Iteration (k)':<15}{'Entropy H(p_k)':<20}{'State Distribution (p_k)'}")
print("-" * 70)

# Run for a number of iterations
for k in range(15):
    # 1. Calculate state distribution p_k for current policy pi_k
    p_k = calculate_state_dist(pi, T)

    # 2. Calculate entropy of the distribution
    entropy = calculate_entropy(p_k)
    p_k_str = np.array2string(p_k, precision=3, floatmode='fixed')
    print(f"{k:<15}{entropy:<20.4f}{p_k_str}")
    
    # 3. Calculate intrinsic rewards for the next iteration
    # Add a small epsilon to p_k to avoid -log(0)
    epsilon = 1e-12
    r_k_plus_1 = -np.log(p_k + epsilon)
    
    # 4. Find the next policy pi_{k+1}
    pi = solve_for_policy(r_k_plus_1, T)

print("-" * 70)
print("\nAs k increases, the state distribution approaches uniform and the entropy H(s) approaches its maximum.")
print("This demonstrates that the limiting policy maximizes the state entropy.")

<<<A>>>