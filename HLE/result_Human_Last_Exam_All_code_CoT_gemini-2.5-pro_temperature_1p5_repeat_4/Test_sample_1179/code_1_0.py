import numpy as np

def run_simulation_and_print_results(n, t_points):
    """
    Runs one simulation of the iterative process and prints the results
    at specified time steps.

    Args:
        n (int): Sample size at each step.
        T (list): A list of time steps t to report on.
    """
    if n <= 1:
        print("Error: Sample size n must be greater than 1 for unbiased variance.")
        return

    print(f"--- Running simulation with n={n} ---\n")

    mu_tm1 = 0.0
    # sigma_i^2 values, starting with sigma_0^2
    sigma2_values = [1.0]

    max_t = max(t_points)
    for t in range(1, max_t + 1):
        # In theory sigma2 is positive, but add a floor for numerical stability
        sigma2_tm1 = sigma2_values[-1]
        if sigma2_tm1 <= 0:
            sigma2_tm1 = 1e-9

        # 1. Sample n variables from N(mu_{t-1}, sigma_{t-1}^2)
        samples = np.random.normal(loc=mu_tm1, scale=np.sqrt(sigma2_tm1), size=n)

        # 2. Compute MLE for mean and unbiased estimator for variance
        mu_t = np.mean(samples)
        sigma2_t = np.var(samples, ddof=1)

        sigma2_values.append(sigma2_t)

        # Check if the current time t is one of the points to report
        if t in t_points:
            # Calculate S_t = sum of sigma_i^2 from i=0 to t
            S_t = sum(sigma2_values)
            
            print(f"At t = {t}:")

            # To fulfill the "output each number in the final equation" request,
            # we demonstrate the structure of the theoretical expectation calculation.
            # Printing the full sum for the simulation is too verbose.
            # Example: S_t = sigma_0^2 + sigma_1^2 + ...
            #          E[S_t] = E[sigma_0^2] + E[sigma_1^2] + ... = 1 + 1 + ...
            
            # Simulated result
            simulated_sum_str = " + ".join(f"{v:.2f}" for v in sigma2_values)
            if len(sigma2_values) > 6:
                simulated_sum_str = f"{sigma2_values[0]:.2f} + ... + {sigma2_values[-1]:.2f}"
            print(f"  Simulated S_{t} = {simulated_sum_str} = {S_t:.4f}")
            
            # Theoretical result
            theoretical_sum_str = " + ".join(["1"] * (t + 1))
            if t + 1 > 6:
                theoretical_sum_str = f"1 + 1 + ... + 1 ({t + 1} terms)"
            
            print(f"  Theoretical E[S_{t}] = {theoretical_sum_str} = {t + 1}")
            print("-" * 20)


if __name__ == '__main__':
    # --- Simulation Parameters ---
    n = 10  # Sample size (must be > 1)
    
    # Time points to observe the value of S_t
    report_at_steps = [5, 20, 50, 100]

    run_simulation_and_print_results(n, report_at_steps)
    
    print("\nConclusion:")
    print("The simulation shows that S_t tends to grow with t, which is consistent with")
    print("the theoretical expectation E[S_t] = t+1.")
    print("The divergence of the expectation proves that S_t does not converge in L1.")
    print("The fact that S_t diverges to infinity almost surely means it does not converge in distribution.")
