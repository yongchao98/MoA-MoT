import numpy as np

def solve_adsorption_model():
    """
    Solves the multi-layer adsorption model using a self-consistent
    mean-field approach based on the provided parameters and assumptions.
    """
    # Parameters from the problem (with assumptions for clarity)
    # Energies are in units of k_B * T
    mu = 0.15
    eps1 = 0.1
    zl = 4
    
    # --- Assumptions for missing/unclear parameters ---
    # 1. Maximum number of layers, k
    k_max = 10
    # 2. Lateral interaction energy, epsilon_l
    # The given formula epsilon_l = (0.02)^k leads to negligible interaction
    # for k > 1. We assume a typo and use the k=1 value.
    epsilon_l = 0.02
    # 3. Upper layer adsorption energy, epsilon_2, ..., epsilon_k
    # We assume epsilon_j for j>1 are all equal to epsilon_2.
    # We estimate epsilon_2 from the cohesive energy, assuming the pair interaction
    # energy is epsilon_l. Total coordination number is (zl + 2).
    epsilon_2 = (zl + 2) / 2.0 * epsilon_l

    print("--- Model Parameters (in units of k_B*T) ---")
    print(f"Chemical Potential (mu): {mu}")
    print(f"First Layer Energy (eps1): {eps1}")
    print(f"Subsequent Layer Energy (eps2): {epsilon_2:.3f} (assumed)")
    print(f"Lateral Interaction (epsilon_l): {epsilon_l} (assumed)")
    print(f"Lateral Coordination (zl): {zl}")
    print(f"Max Layers (k_max): {k_max} (assumed)")
    print("-------------------------------------------\n")


    # Pack layer energies into an array for convenience (use 1-based indexing)
    epsilon_j = np.zeros(k_max + 1)
    epsilon_j[1] = eps1
    epsilon_j[2:] = epsilon_2

    # --- Self-consistent iterative solver for theta_j ---
    # Initial guess: no particles adsorbed
    thetas = np.zeros(k_max + 1)
    
    # Iteration control
    max_iter = 1000
    tolerance = 1e-7
    mixing_alpha = 0.1  # Mixing factor to ensure stable convergence

    for i in range(max_iter):
        thetas_old = np.copy(thetas)
        
        # 1. Calculate effective energies based on current thetas
        epsilon_eff = np.zeros(k_max + 1)
        epsilon_eff[1] = epsilon_j[1] + zl * epsilon_l * thetas[1]
        for j in range(2, k_max + 1):
            epsilon_eff[j] = epsilon_j[j] + zl * epsilon_l * thetas[j]

        # 2. Calculate terms for the partition function Xi = 1 + sum(x_m)
        # where x_m = exp(m*mu + sum_{j=1 to m} epsilon_eff_j)
        # To avoid numerical overflow, we work with logarithms.
        log_x = np.zeros(k_max + 1)
        current_sum_eps = 0
        for m in range(1, k_max + 1):
            current_sum_eps += epsilon_eff[m]
            log_x[m] = m * mu + current_sum_eps
        
        # Normalize to prevent overflow: find max term, divide by it, then multiply back.
        log_max_val = np.max(np.concatenate(([0], log_x[1:])))
        
        # Calculate terms relative to the max term
        xi_terms_norm = np.exp(log_x[1:] - log_max_val)
        xi_site_norm = np.exp(0 - log_max_val) + np.sum(xi_terms_norm)
        
        # 3. Calculate new thetas
        # theta_j = Sum_{m=j to k} P(m) = (Sum_{m=j to k} x_m) / Xi
        thetas_new = np.zeros(k_max + 1)
        cumulative_sum_terms = 0
        # Iterate backwards from j=k_max to efficiently calculate the sum
        for j in range(k_max, 0, -1):
            cumulative_sum_terms += xi_terms_norm[j - 1] # 0-indexed array
            thetas_new[j] = cumulative_sum_terms / xi_site_norm
        
        # 4. Mix new and old thetas for stability
        thetas = mixing_alpha * thetas_new + (1 - mixing_alpha) * thetas_old

        # 5. Check for convergence
        if np.max(np.abs(thetas - thetas_old)) < tolerance:
            print(f"Converged after {i+1} iterations.\n")
            break
    else:
        print("Warning: Did not converge within max iterations.\n")

    # --- Final Calculation ---
    # The average number of layers <k> is the sum of the layer coverages.
    avg_k = np.sum(thetas[1:])

    print("Final converged layer coverages:")
    theta_strings = []
    for j in range(1, k_max + 1):
        print(f"  theta_{j:2d}: {thetas[j]:.6f}")
        theta_strings.append(f"{thetas[j]:.6f}")
    
    print("\nThe average number of adsorbed layers per site, <k>, is calculated as:")
    print(f"<k> = sum(theta_j) for j=1..{k_max}")
    print(f"<k> = {' + '.join(theta_strings)}")
    print(f"\nFinal Result: <k> = {avg_k:.6f}")
    return avg_k

if __name__ == '__main__':
    final_avg_k = solve_adsorption_model()
    # The final answer is wrapped according to the required format
    print(f"\n<<<{final_avg_k:.4f}>>>")
