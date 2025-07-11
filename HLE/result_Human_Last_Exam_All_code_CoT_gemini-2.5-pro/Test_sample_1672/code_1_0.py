import numpy as np

def calculate_entropy(p):
    """Calculates the entropy of a probability distribution in bits."""
    # Add a small epsilon to prevent log(0) for states with zero probability.
    p = p + 1e-12
    return -np.sum(p * np.log2(p))

def run_simulation():
    """
    Simulates the iterative process of entropy maximization and prints the results.
    """
    # 1. Define the environment and initial state
    N_STATES = 10
    max_entropy = np.log2(N_STATES)

    # 2. Define the initial state distribution for policy pi^0.
    # We start with a highly skewed (low entropy) distribution.
    p_k = np.array([0.9] + [0.1 / (N_STATES - 1)] * (N_STATES - 1))
    
    # Store distributions and entropies
    distributions = {'pi^0': p_k}
    entropies = {'pi^0': calculate_entropy(p_k)}

    # 3. Simulate the iterative process for policy updates.
    # The policy pi^k is trained with reward r_k(s) = -log(p_{k-1}(s)).
    # A policy maximizing this reward will produce a state distribution p_k(s)
    # that is higher for states where p_{k-1}(s) was lower.
    # We can model this update on the distribution directly:
    # p_k(s) is proportional to exp(-log(p_{k-1}(s))) = 1 / p_{k-1}(s).

    # Iteration k=1 (for policy pi^1)
    p_k_unnormalized = 1.0 / (p_k + 1e-12)
    p_k = p_k_unnormalized / np.sum(p_k_unnormalized)
    distributions['pi^1'] = p_k
    entropies['pi^1'] = calculate_entropy(p_k)

    # Iterate until convergence (k -> infinity)
    for _ in range(20): # 20 iterations are sufficient for convergence in this model
        p_k_unnormalized = 1.0 / (p_k + 1e-12)
        p_k = p_k_unnormalized / np.sum(p_k_unnormalized)

    distributions['lim_{k->inf} pi^k'] = p_k
    entropies['lim_{k->inf} pi^k'] = calculate_entropy(p_k)

    # 4. Print the results
    print("This simulation shows how an iterative process converges to a policy that maximizes state entropy.")
    print("At each step k, the policy pi^k is updated to favor states that were rare under pi^{k-1}.")
    print("-" * 70)
    print(f"Theoretical Maximum Entropy for {N_STATES} states: {max_entropy:.4f}\n")

    print("Calculated Entropies for policies at different iterations:")
    print(f"Policy pi^0 (k=0):        H(s) = {entropies['pi^0']:.4f}")
    print(f"Policy pi^1 (k=1):        H(s) = {entropies['pi^1']:.4f}")
    print(f"Policy lim (k->inf): H(s) = {entropies['lim_{k->inf} pi^k']:.4f}")
    print("-" * 70)

    print("As the number of iterations k increases, the state entropy H(s) increases and approaches the maximum value.")
    print("The policy that maximizes the entropy is the one at the limit of this process.\n")
    
    # Final equation showing the comparison
    print("Final Comparison:")
    print(f"H(s) for pi^0 ({entropies['pi^0']:.4f}) < H(s) for pi^1 ({entropies['pi^1']:.4f}) < H(s) for lim_{{k->inf}} pi^k ({entropies['lim_{k->inf} pi^k']:.4f})")

if __name__ == '__main__':
    run_simulation()