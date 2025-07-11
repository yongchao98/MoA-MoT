import math

def calculate_lowest_coverage(n, alpha):
    """
    Calculates the guaranteed lowest coverage probability for a LOO-based
    conformal prediction interval.

    Args:
        n (int): The number of training data points.
        alpha (float): The significance level (e.g., 0.1 for 90% desired coverage).

    Returns:
        float: The minimum guaranteed coverage probability.
    """
    if not isinstance(n, int) or n <= 0:
        raise ValueError("n must be a positive integer.")
    if not isinstance(alpha, float) or not (0 < alpha < 1):
        raise ValueError("alpha must be a float between 0 and 1.")

    # The number of scores used for the quantile calculation is n+1
    # (n LOO scores + infinity)
    n_plus_1 = n + 1

    # The quantile is based on the level 1-alpha
    one_minus_alpha = 1 - alpha

    # The rank k is the ceiling of (n+1)*(1-alpha)
    k_numerator = n_plus_1 * one_minus_alpha
    k = math.ceil(k_numerator)

    # The lowest guaranteed coverage probability is k / (n+1)
    coverage = k / n_plus_1

    # --- Printing the detailed calculation ---
    print(f"For n = {n} and alpha = {alpha}:")
    print("The lowest guaranteed coverage probability is given by the formula:")
    print("P(Coverage) = ceil((n + 1) * (1 - alpha)) / (n + 1)\n")
    print("Step-by-step calculation:")
    print(f"1. n + 1 = {n} + 1 = {n_plus_1}")
    print(f"2. 1 - alpha = 1 - {alpha} = {one_minus_alpha}")
    print(f"3. (n + 1) * (1 - alpha) = {n_plus_1} * {one_minus_alpha} = {k_numerator}")
    print(f"4. k = ceil({k_numerator}) = {k}")
    print(f"5. Lowest Coverage = k / (n + 1) = {k} / {n_plus_1} = {coverage}")
    print("-" * 20)

    return coverage

# --- Example Usage ---
# Example 1: n=19, alpha=0.1. (1-alpha)*(n+1) is an integer.
calculate_lowest_coverage(n=19, alpha=0.1)

# Example 2: n=20, alpha=0.1. (1-alpha)*(n+1) is not an integer.
calculate_lowest_coverage(n=20, alpha=0.1)

# Example 3: A more typical case with n=999 and alpha=0.05
calculate_lowest_coverage(n=999, alpha=0.05)
