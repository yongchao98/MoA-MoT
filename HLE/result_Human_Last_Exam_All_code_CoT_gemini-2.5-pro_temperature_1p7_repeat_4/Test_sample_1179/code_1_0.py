import numpy as np

def simulate_variance_sum(n, num_steps):
    """
    Simulates the iterative process of sampling from a normal distribution
    and updating its parameters, while tracking the sum of variances.

    Args:
        n (int): The number of samples to draw at each step.
        num_steps (int): The total number of steps (t) to run the simulation.
    """
    if n < 2:
        print("Error: n must be 2 or greater for the unbiased variance to be defined.")
        return

    # Initial parameters at t=0
    mu = 0.0
    sigma2 = 1.0

    # Initial sum S_0
    S = sigma2
    
    print(f"Simulation with n = {n}\n")
    # Print the state at t=0
    print(f"S_0 = sigma^2_0 = {sigma2:.4f}")

    previous_S = S
    for t in range(1, num_steps + 1):
        # Ensure variance is non-negative before taking the square root for scale
        scale = np.sqrt(max(0, sigma2))
        
        # Step 1: Sample n variables from N(mu_{t-1}, sigma^2_{t-1})
        samples = np.random.normal(loc=mu, scale=scale, size=n)

        # Step 2: Compute new estimators
        # mu_t is the sample mean
        mu_t = np.mean(samples)
        # sigma^2_t is the unbiased sample variance (ddof=1)
        sigma2_t = np.var(samples, ddof=1)
        
        # Update the sum S_t = S_{t-1} + sigma^2_t
        S += sigma2_t

        # Print the equation for the current step's sum
        print(f"S_{t} = S_{t-1} + sigma^2_{t} = {previous_S:.4f} + {sigma2_t:.4f} = {S:.4f}")

        # Update parameters for the next iteration
        mu = mu_t
        sigma2 = sigma2_t
        previous_S = S

# --- Simulation Parameters ---
# Number of samples per step
n = 5
# Total number of iterations
T = 25

# Run the simulation
simulate_variance_sum(n, T)