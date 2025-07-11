import numpy as np

def run_simulation():
    """
    Runs a simulation to demonstrate the behavior of S_t.

    The function simulates the described iterative process multiple times to
    numerically verify the theoretical expectation of S_t.
    """
    # --- Parameters ---
    n = 5  # Sample size at each step (must be >= 2)
    num_simulations = 20000  # Number of independent trials to run
    
    # Time points at which we will check the value of S_t
    checkpoints = [10, 50, 100, 200]
    max_t = max(checkpoints)

    # Dictionary to store the list of S values at each checkpoint
    s_values_at_checkpoints = {t: [] for t in checkpoints}

    print(f"Running {num_simulations} simulations with n={n}...")
    
    # --- Simulation Loop ---
    for i in range(num_simulations):
        # Initial conditions for each simulation trial
        mu_t = 0.0
        sigma2_t = 1.0
        
        # S_t starts with S_0 = sigma_0^2
        s_t = sigma2_t

        # Loop through time steps
        for t in range(1, max_t + 1):
            # 1. Sample n variables from N(mu_{t-1}, sigma^2_{t-1})
            # Handle cases where sigma2_t becomes very small or zero
            if sigma2_t <= 0:
                # This is unlikely but can happen due to floating point inaccuracies.
                # If it happens, we stop this trial's evolution.
                break 
            samples = np.random.normal(loc=mu_t, scale=np.sqrt(sigma2_t), size=n)
            
            # 2. Compute new mu_t and sigma2_t
            mu_t = np.mean(samples)
            # Use ddof=1 for the unbiased estimator of the variance
            sigma2_t = np.var(samples, ddof=1)
            
            # 3. Update the sum S_t
            s_t += sigma2_t
            
            # 4. Store the result if it's a checkpoint
            if t in checkpoints:
                s_values_at_checkpoints[t].append(s_t)

    print("\n--- Results ---")
    print("Comparing simulated average of S_t with theoretical expectation E[S_t] = t + 1.\n")

    # --- Output Results ---
    for t in checkpoints:
        simulated_mean_s = np.mean(s_values_at_checkpoints[t])
        theoretical_mean_s = float(t + 1)
        
        print(f"Time step t = {t}:")
        print(f"  Theoretical Expectation E[S_t] = {t} + 1 = {theoretical_mean_s}")
        print(f"  Simulated Average S_t          = {simulated_mean_s:.4f}")
        print("-" * 20)

if __name__ == '__main__':
    run_simulation()
