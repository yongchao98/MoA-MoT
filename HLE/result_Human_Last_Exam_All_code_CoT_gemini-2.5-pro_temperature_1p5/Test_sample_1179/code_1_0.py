import numpy as np

def simulate_convergence(n, num_paths, max_t):
    """
    Simulates the iterative process to observe the convergence of S_t.

    Args:
        n (int): The number of samples drawn at each step.
        num_paths (int): The number of independent simulations to run.
        max_t (int): The maximum number of steps in the simulation.
    """
    if n <= 1:
        print("Error: n must be greater than 1 for the unbiased variance to be defined.")
        return

    # Initial conditions for all paths
    # sigma2_t stores the variance sigma_t^2 for each path
    # We start with t=0, so sigma_0^2 = 1
    sigma2_t = np.ones(num_paths)
    
    # S_t stores the sum sum_{i=0 to t} sigma_i^2 for each path
    # We start with S_0 = sigma_0^2 = 1
    S_t = np.ones(num_paths)

    # Degrees of freedom for the chi-squared distribution
    df = n - 1

    print(f"Simulation parameters: n={n}, num_paths={num_paths}")
    print("-" * 70)
    print(f"{'Time (t)':>10s} | {'Mean(S_t)':>15s} | {'Std Dev(S_t)':>15s} | {'P50(S_t)':>10s} | {'P95(S_t)':>10s}")
    print("-" * 70)

    checkpoints = [0, 10, 50, 100, 200, 500, 1000, 2000, 5000]

    for t in range(1, max_t + 1):
        # Generate W_t = chi^2(n-1) / (n-1) for all paths
        W_t = np.random.chisquare(df, num_paths) / df
        
        # Update sigma_t^2 = sigma_{t-1}^2 * W_t
        sigma2_t *= W_t
        
        # Update S_t = S_{t-1} + sigma_t^2
        S_t += sigma2_t

        # Print statistics at specified checkpoints
        if t in checkpoints:
            mean_S = np.mean(S_t)
            std_S = np.std(S_t)
            p50_S = np.percentile(S_t, 50) # Median
            p95_S = np.percentile(S_t, 95)
            print(f"{t:>10d} | {mean_S:>15.2f} | {std_S:>15.2f} | {p50_S:>10.2f} | {p95_S:>10.2f}")

if __name__ == '__main__':
    # --- Parameters ---
    # Sample size at each step
    n_samples = 4
    # Number of simulations to get a distribution
    num_simulations = 50000
    # Maximum number of steps
    max_steps = 5000
    
    simulate_convergence(n_samples, num_simulations, max_steps)