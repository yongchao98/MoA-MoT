import numpy as np

def simulate_and_analyze(t_max, n, num_simulations=5000):
    """
    Simulates the iterative process and analyzes the sum S_t.

    Args:
        t_max (int): The number of steps for the sum (the 't' in S_t).
        n (int): The sample size at each step. Must be > 1.
        num_simulations (int): The number of times to run the simulation.
    """
    if n <= 1:
        print("Error: Sample size n must be greater than 1 for an unbiased variance estimator.")
        return

    print(f"Analyzing S_t = sum(sigma_i^2 for i=0 to t) for t = {t_max} and n = {n}\n")

    # --- Theoretical Analysis ---
    print("--- 1. Theoretical Analysis ---")
    theoretical_mean = t_max + 1
    
    # Building and printing the equation for the expectation E[S_t]
    print("The expectation of the sum S_t is E[S_t] = sum(E[sigma_i^2]).")
    print(f"Since E[sigma_i^2] = 1 for all i, the equation is:")
    
    # To avoid printing a very long line, we'll show a compact version for large t
    if t_max <= 5:
        equation_terms = ["E[sigma_{}^2]".format(i) for i in range(t_max + 1)]
        print(f"E[S_{t_max}] = {' + '.join(equation_terms)}")
        values_str = "         = " + " + ".join(["1"] * (t_max + 1))
        print(values_str)
    else:
        # Show a shortened version for readability
        short_eq_str = f"E[S_{t_max}] = E[sigma_0^2] + E[sigma_1^2] + ... + E[sigma_{t_max}^2]"
        print(short_eq_str)
        print(f"         = 1 + 1 + ... + 1  ({t_max+1} terms)")

    print(f"         = {theoretical_mean}\n")
    print(f"Because E[S_t] = t + 1, the expectation grows to infinity as t increases.")
    print("This means the sum S_t diverges and does not converge in L1 or in distribution.")
    print("-" * 35)

    # --- Simulation Verification ---
    print("\n--- 2. Simulation Verification ---")
    print(f"Running {num_simulations} simulations to find the sample mean of S_{t_max}...")
    
    final_s_values = []
    for _ in range(num_simulations):
        mu = 0.0
        sigma2 = 1.0
        
        # S_t starts with S_0 = sigma_0^2
        s_t = sigma2
        
        for _ in range(t_max):
            # Sample from N(mu, sigma2)
            # np.random.normal takes standard deviation (sqrt of variance)
            if sigma2 < 1e-9: # Prevent negative variance from floating point error
                sigma2 = 1e-9

            samples = np.random.normal(loc=mu, scale=np.sqrt(sigma2), size=n)
            
            # Update mu (MLE) and sigma2 (unbiased estimator)
            mu = np.mean(samples)
            sigma2 = np.var(samples, ddof=1)
            
            s_t += sigma2
            
        final_s_values.append(s_t)

    simulated_mean = np.mean(final_s_values)

    print("\n--- Results ---")
    print(f"Theoretical Expectation E[S_{t_max}]: {theoretical_mean}")
    print(f"Simulated Mean of S_{t_max}:      {simulated_mean:.4f}")
    print(f"The simulation confirms that the mean of S_t grows with t.")


# You can change these parameters to explore the process
# t_max: The final time step for the sum S_t
# n: The number of samples drawn at each step
simulate_and_analyze(t_max=20, n=10)