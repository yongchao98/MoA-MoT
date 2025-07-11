import numpy as np

def run_simulation():
    """
    Simulates the iterative process to analyze the convergence of S_t.
    """
    # --- Simulation Parameters ---
    n = 10                # Sample size at each step
    num_steps = 101       # Total number of steps (from t=0 to t=100)
    num_trajectories = 5000 # Number of independent simulations to run

    # --- Initial Conditions ---
    mu_0 = 0.0
    sigma2_0 = 1.0

    # Array to store the history of the sum S_t for each trajectory
    s_history = np.zeros((num_trajectories, num_steps))

    print(f"Running {num_trajectories} simulations for {num_steps-1} steps with n={n}...\n")

    # Loop over each independent trajectory
    for i in range(num_trajectories):
        mu = mu_0
        sigma2 = sigma2_0
        
        current_s = sigma2_0
        s_history[i, 0] = current_s

        # Loop over each time step t > 0
        for t in range(1, num_steps):
            # Ensure variance is non-negative before sampling
            if sigma2 < 0:
                sigma2 = 0
            
            # 1. Sample n variables from N(mu_{t-1}, sigma2_{t-1})
            samples = np.random.normal(loc=mu, scale=np.sqrt(sigma2), size=n)

            # 2. Compute new mu_t and sigma2_t
            # MLE for mean is the sample mean
            mu_t = np.mean(samples)
            # Unbiased estimator for variance is the sample variance (ddof=1)
            sigma2_t = np.var(samples, ddof=1)
            
            # Handle rare case where sample variance is 0
            if np.isnan(sigma2_t):
                sigma2_t = 0
            
            # Update parameters for the next step
            mu = mu_t
            sigma2 = sigma2_t

            # Update the sum S_t
            current_s += sigma2
            s_history[i, t] = current_s

    # --- Analysis ---
    # Calculate the mean of S_t across all trajectories at each step
    mean_s_t = np.mean(s_history, axis=0)

    print("Comparing empirical average of S_t with theoretical expectation (t+1):")
    print("-" * 70)
    print(f"{'Time Step (t)':<15} | {'Empirical E[S_t]':<25} | {'Theoretical E[S_t] = t+1':<25}")
    print("-" * 70)

    # Print results at intervals
    for t in range(0, num_steps, 10):
        empirical_mean = mean_s_t[t]
        theoretical_mean = t + 1
        print(f"{t:<15} | {empirical_mean:<25.4f} | {theoretical_mean:<25}")

    print("-" * 70)
    print("\nAs the simulation shows, the average value of S_t grows linearly with t,")
    print("supporting the theoretical conclusion that the sum diverges.")


if __name__ == '__main__':
    run_simulation()