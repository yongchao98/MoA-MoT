import numpy as np

def simulate_and_analyze():
    """
    Simulates the described iterative process and analyzes the convergence of S_t.
    """
    # Parameters for the simulation
    n = 10                # Sample size at each step
    T = 100               # Total number of iterations (time steps)
    num_simulations = 5000  # Number of full simulations to run for averaging

    # Define time points to check the results
    checkpoints = [10, 50, T]
    
    # Store the final S_t value for each simulation at each checkpoint
    results = {t: [] for t in checkpoints}

    # Run the simulation multiple times
    for _ in range(num_simulations):
        # Initialization for each simulation run (t=0)
        mu = 0.0
        sigma2 = 1.0
        S_t = sigma2 # S_0 = sigma_0^2 = 1

        # Iterate from t=1 to T
        for t in range(1, T + 1):
            # 1. Sample n variables from N(mu_{t-1}, sigma_{t-1}^2)
            # A small floor for sigma2 to prevent numerical issues if it gets too close to zero
            if sigma2 < 1e-12:
                sigma2 = 1e-12
            std_dev = np.sqrt(sigma2)
            samples = np.random.normal(loc=mu, scale=std_dev, size=n)
            
            # 2. Compute new mu_t and sigma_t^2
            mu = np.mean(samples)
            # Use ddof=1 for the unbiased estimator of the variance
            sigma2 = np.var(samples, ddof=1)
            
            # 3. Update the sum S_t
            S_t += sigma2
            
            # Store the sum if the current time t is a checkpoint
            if t in checkpoints:
                results[t].append(S_t)

    # Print the analysis
    print("--- Analysis of S_t = sum(sigma_i^2) for i=0 to t ---")
    print("\nTheoretical expectation E[S_t] is t + 1.")
    print("Let's compare this with the simulated mean of S_t.\n")
    
    for t in checkpoints:
        theoretical_mean_eq = f"{t} + 1 = {t+1}"
        simulated_mean = np.mean(results[t])
        
        print(f"Time t = {t}:")
        print(f"  - Theoretical Expectation: E[S_{t}] = {theoretical_mean_eq}")
        print(f"  - Simulated Mean of S_{t}: {simulated_mean:.4f}")

    print("\n--- Conclusion ---")
    print("The expectation E[S_t] grows linearly with t and goes to infinity.")
    print("This means S_t cannot converge in L1.")
    print("\nFurthermore, the sequence of variances, sigma_t^2, forms a martingale that converges")
    print("almost surely to a random variable with mean 1, not 0. Since the terms of the")
    print("series S_t do not approach 0, the series diverges almost surely to infinity.")
    print("A sequence that diverges to infinity does not converge in distribution to a proper random variable.")

simulate_and_analyze()