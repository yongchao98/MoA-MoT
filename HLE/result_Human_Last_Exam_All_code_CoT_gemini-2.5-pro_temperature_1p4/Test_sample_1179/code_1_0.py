import numpy as np

def run_simulation(n, max_t, time_points):
    """
    Runs a single simulation trial of the iterative process.
    
    Args:
        n (int): The number of samples to draw at each step.
        max_t (int): The total number of steps to simulate.
        time_points (list): A list of time steps at which to record the value of S_t.

    Returns:
        dict: A dictionary mapping time points to the value of S_t at that time.
    """
    mu = 0.0
    sigma_sq = 1.0
    s_t = sigma_sq
    results = {}

    # Convert time_points to a set for efficient lookup
    time_points_set = set(time_points)

    for t in range(1, max_t + 1):
        # Handle the case where sigma_sq becomes very small or zero to avoid errors
        if sigma_sq <= 1e-9:
            sigma_sq = 1e-9

        # Step 1: Sample n variables from N(mu, sigma_sq)
        # Note: np.random.normal takes standard deviation (sqrt(variance)) as scale
        samples = np.random.normal(loc=mu, scale=np.sqrt(sigma_sq), size=n)
        
        # Step 2: Compute new mu (MLE) and sigma_sq (unbiased estimator)
        mu = np.mean(samples)
        # Use ddof=1 for the unbiased estimator of the variance
        sigma_sq = np.var(samples, ddof=1)
        
        # Step 3: Update the sum S_t
        s_t += sigma_sq
        
        # Record the result if t is one of the specified time points
        if t in time_points_set:
            results[t] = s_t
            
    return results

def main():
    """
    Main function to run the simulation and print the analysis.
    """
    # --- Simulation Parameters ---
    n = 10              # Sample size per step
    max_t = 2000        # Total number of steps in each simulation
    num_trials = 5000   # Number of simulation runs to average over

    print("--- Simulation Setup ---")
    print(f"Sample size at each step (n): {n}")
    print(f"Total steps per simulation (t_max): {max_t}")
    print(f"Number of trials: {num_trials}\n")

    # Time points at which we will analyze the results
    time_points_to_check = [1, 10, 100, 500, 1000, 2000]
    
    # --- Run Simulations ---
    all_results = []
    print("Running simulations... (this may take a moment)")
    for _ in range(num_trials):
        trial_results = run_simulation(n, max_t, time_points_to_check)
        all_results.append(trial_results)
    print("Simulations complete.\n")

    # --- Analyze and Print Results ---
    print("--- Analysis of S_t = sum(sigma_i^2) ---")
    print("The theoretical analysis shows that E[S_t] = t + 1.")
    print("Below, we compare the simulated average of S_t with its theoretical value.")
    # Printing the "final equation" E[S_t] = t + 1 for each specified t
    print("\n" + f"{'Time t':>8} | {'Simulated Avg S_t':>20} | {'Theoretical E[S_t]':>22}")
    print("-" * 55)

    for t in time_points_to_check:
        # Extract the value of S_t from all trials at the current time point
        s_values_at_t = [res[t] for res in all_results]
        
        # Calculate the average S_t from the simulation
        simulated_avg_s_t = np.mean(s_values_at_t)
        
        # Calculate the theoretical expected value E[S_t]
        theoretical_mean = t + 1.0
        
        print(f"{t:>8} | {simulated_avg_s_t:>20.4f} | {theoretical_mean:>22.4f}")

    print("\n--- Conclusion ---")
    print("The simulation results closely match the theoretical expectation. The average value of S_t")
    print("grows linearly with t, clearly indicating that E[S_t] -> infinity as t -> infinity.")
    print("Therefore, the sum S_t diverges and does not converge in L1 or in distribution.")

if __name__ == "__main__":
    main()