import numpy as np

def simulate_variance_sum(n, num_steps):
    """
    Simulates the iterative process for variance and calculates the sum S_t.
    This function demonstrates that the sum S_t grows over time and does not converge.

    The final "equation" for S_t is shown at each step by printing the
    newly calculated sigma_t^2 and the resulting sum S_t.

    Args:
        n (int): The number of samples at each step (must be >= 2).
        num_steps (int): The total number of steps (t) to simulate.
    """
    if n < 2:
        print("Error: Sample size n must be 2 or greater.")
        return

    # Initial conditions at t=0
    mu_prev = 0.0
    sigma2_prev = 1.0

    # Initialize S_t with the t=0 term: S_0 = sigma_0^2
    s_t = sigma2_prev
    
    # The output shows each number contributing to the final sum
    print(f"Initial state (t=0):")
    print(f"sigma_0^2 = {sigma2_prev:.4f}")
    print(f"S_0 = {s_t:.4f}")
    print("-" * 30)

    # Iterative process for t = 1, 2, ..., num_steps
    for t in range(1, num_steps + 1):
        # The scale parameter for numpy's normal distribution is the standard deviation
        std_dev_prev = np.sqrt(sigma2_prev)
        
        # Sample n variables from N(mu_{t-1}, sigma_{t-1}^2)
        samples = np.random.normal(loc=mu_prev, scale=std_dev_prev, size=n)

        # Compute mu_t (MLE for mean) and sigma_t^2 (unbiased estimator for variance)
        mu_t = np.mean(samples)
        # Using ddof=1 provides the unbiased sample variance
        sigma2_t = np.var(samples, ddof=1)
        
        # In rare cases of floating-point error, variance could be slightly negative.
        if sigma2_t < 0:
            sigma2_t = 0

        # Update the sum S_t = S_{t-1} + sigma_t^2
        s_t += sigma2_t
        
        print(f"Step t={t}:")
        print(f"Sampled from N(mean={mu_prev:.4f}, var={sigma2_prev:.4f})")
        # Output the numbers for the "equation" S_t = S_{t-1} + sigma_t^2
        print(f"sigma_{t}^2 = {sigma2_t:.4f}")
        print(f"S_{t} = S_{t-1} + sigma_{t}^2 = {s_t:.4f}")
        print("-" * 30)

        # Update values for the next iteration
        mu_prev = mu_t
        sigma2_prev = sigma2_t

# --- Parameters for the simulation ---
# Sample size at each step (n). Must be >= 2.
sample_size = 10
# Total number of steps to simulate (t_max).
total_steps = 25

# Run the simulation
simulate_variance_sum(sample_size, total_steps)