import numpy as np

def simulate_process():
    """
    This function simulates the iterative process and demonstrates the divergence 
    of the sum S_t by comparing its simulated average with the theoretical expectation.
    """
    # --- Parameters ---
    n = 10                # Sample size at each step
    max_t = 100           # Maximum number of iterations
    num_simulations = 20000 # Number of simulations to average over

    # Array to store the sum of S_t over all simulations for each t
    # s_t_accumulator[t] will store the sum of S_t values from all simulations
    s_t_accumulator = np.zeros(max_t + 1)

    print(f"Running {num_simulations} simulations up to t={max_t} with n={n}...\n")

    # --- Simulation Loop ---
    for _ in range(num_simulations):
        # Initialize the process for a single simulation run
        mu = 0.0
        sigma2 = 1.0
        
        # S_t for the current simulation run, starting with S_0 = sigma_0^2
        current_s_t = sigma2
        
        # Accumulate S_0
        s_t_accumulator[0] += current_s_t

        for t in range(1, max_t + 1):
            # Step 1: Sample n variables from N(mu, sqrt(sigma2))
            # np.random.normal requires standard deviation (sqrt of variance)
            # Handle potential floating point inaccuracies where sigma2 might be slightly negative
            std_dev = np.sqrt(max(0, sigma2))
            samples = np.random.normal(loc=mu, scale=std_dev, size=n)

            # Step 2: Compute new mu_t (MLE) and sigma_t^2 (unbiased)
            mu = np.mean(samples)
            # Use ddof=1 for the unbiased estimator of the variance
            sigma2 = np.var(samples, ddof=1)
            
            # Step 3: Update the sum S_t
            current_s_t += sigma2
            
            # Accumulate the S_t value for the current t
            s_t_accumulator[t] += current_s_t

    # --- Analysis and Output ---
    # Calculate the average S_t across all simulations
    avg_s_t = s_t_accumulator / num_simulations

    print("Comparing theoretical E[S_t] with simulated average E[S_t]:")
    print("-" * 65)
    print(f"{'Time t':<10} | {'Theoretical E[S_t] = t + 1':<30} | {'Simulated E[S_t]':<20}")
    print("-" * 65)

    # Time points to display
    time_points = [0, 10, 20, 50, 80, max_t]

    for t in time_points:
        theoretical_exp = t + 1
        simulated_exp = avg_s_t[t]
        
        # Print the equation with its numbers
        print(f"t = {t:<7} | E[S_{t}] = {t} + 1 = {theoretical_exp:<22.2f} | {simulated_exp:<20.2f}")
        
    print("-" * 65)
    print("\nAs shown, the expected value of S_t grows linearly with t and does not converge.")
    print("This supports the theoretical conclusion that S_t diverges.")

# Run the simulation
simulate_process()
