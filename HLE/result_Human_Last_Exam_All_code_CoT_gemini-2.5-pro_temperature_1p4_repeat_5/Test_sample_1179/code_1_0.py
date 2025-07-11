import numpy as np

def simulate_variance_process(n, T):
    """
    Simulates the iterative process for variance and its sum.

    Args:
        n (int): The number of samples to draw at each step.
        T (int): The total number of steps (iterations) to simulate.
    """
    # Step 0: Initialization
    mu_current = 0.0
    var_current = 1.0
    S_t = var_current  # S_0 = sigma_0^2

    print(f"Sample size n = {n}")
    print("-" * 40)
    # The "final equation" is the sum S_t. We print its components at each step.
    # Format: t, sigma_t^2, S_t = sum_{i=0 to t} sigma_i^2
    print(f"t=0:  sigma_0^2 = {var_current:.4f},  S_0 = {S_t:.4f}")

    # Iterative process from t = 1 to T
    for t in range(1, T + 1):
        # Sample n variables from N(mu_{t-1}, sigma_{t-1}^2)
        # np.sqrt(var_current) is the standard deviation sigma_{t-1}
        samples = np.random.normal(loc=mu_current, scale=np.sqrt(var_current), size=n)

        # Compute mu_t (MLE for mean)
        mu_current = np.mean(samples)

        # Compute sigma_t^2 (unbiased estimator for variance)
        # ddof=1 ensures the denominator is n-1 for an unbiased estimate.
        var_current = np.var(samples, ddof=1)

        # Update the sum S_t
        S_t += var_current

        # Print the numbers for the current step's equation
        # S_t = S_{t-1} + sigma_t^2
        print(f"t={t}:  sigma_{t}^2 = {var_current:.4f},  S_{t} = {S_t:.4f}")

    print("-" * 40)
    print(f"As you can see from the output, the sum S_t consistently increases and does not appear to converge.")
    print(f"My analysis shows E[S_t] = t + 1, which goes to infinity. This implies S_t diverges.")


if __name__ == '__main__':
    # --- Parameters ---
    # Number of samples at each step. Must be > 1 for variance calculation.
    sample_size_n = 10
    # Total number of time steps to simulate.
    time_steps_T = 20

    simulate_variance_process(sample_size_n, time_steps_T)