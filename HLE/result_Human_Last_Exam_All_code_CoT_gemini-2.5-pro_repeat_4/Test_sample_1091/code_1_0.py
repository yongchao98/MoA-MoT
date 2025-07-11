import numpy as np

def estimate_limit(k_val, num_trials):
    """
    Numerically estimates the limit of n*P(n) using Monte Carlo simulation.

    Args:
        k_val (int): A large integer to define the number of vectors.
        num_trials (int): The number of simulations to run.
    """
    # n is the total number of vectors, m is the number of vectors of each type.
    n = 6 * k_val
    m = 2 * k_val

    # For efficiency, we can approximate the sum of Rademacher variables
    # with a normal distribution, as per the Central Limit Theorem.
    # X ~ Sum of m Rademacher variables is approx N(0, m).
    # We generate num_trials samples for X_A, X_B, and X_C at once.
    std_dev = np.sqrt(m)
    X_A = np.random.normal(0, std_dev, num_trials)
    X_B = np.random.normal(0, std_dev, num_trials)
    X_C = np.random.normal(0, std_dev, num_trials)

    # Calculate the components of the sum vector S for all trials
    S_x = X_A + 0.5 * (X_B - X_C)
    S_y = (np.sqrt(3) / 2) * (X_B + X_C)

    # Calculate the squared norm of S for all trials
    S_norm_sq = S_x**2 + S_y**2

    # Count the number of "successes" where the norm condition is met
    # The condition is ||S||_2 <= sqrt(2), which is ||S||_2^2 <= 2
    success_count = np.sum(S_norm_sq <= 2)

    # Estimate the probability P(n)
    P_n = success_count / num_trials

    # Calculate the final value n * P(n)
    limit_estimate = n * P_n

    print(f"Simulation parameters:")
    print(f"k = {k_val}, n = {n}, m = {m}")
    print(f"Number of trials = {num_trials}")
    print("-" * 20)
    print(f"Number of successes (||S||^2 <= 2): {success_count}")
    print(f"Estimated P(n) = {P_n:.4g}")
    print(f"Final estimated value for n * P(n) is: {limit_estimate:.4f}")
    print("\nThe analytical result for the limit is 2.")

if __name__ == '__main__':
    # We choose a large k to approximate the limit n -> infinity
    # and a large number of trials for an accurate estimate.
    k_value = 20000
    number_of_trials = 1000000
    estimate_limit(k_value, number_of_trials)
