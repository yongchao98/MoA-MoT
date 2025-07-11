import math

def calculate_lowest_coverage(n, alpha):
    """
    Calculates the theoretical lower bound for the coverage probability
    of a leave-one-out based conformal prediction interval.

    Args:
        n (int): The number of training data points.
        alpha (float): The desired miscoverage level (e.g., 0.1 for 90% confidence).

    Returns:
        float: The guaranteed lowest coverage probability.
    """
    if not isinstance(n, int) or n <= 0:
        raise ValueError("n must be a positive integer.")
    if not (0 < alpha < 1):
        raise ValueError("alpha must be between 0 and 1.")

    print(f"Calculating for n = {n} and alpha = {alpha}")

    # The total number of observations when including the test point
    N = n + 1
    print(f"Total number of points (n+1): {N}")

    # The target coverage level
    target_quantile_level = 1 - alpha
    print(f"Target quantile level (1-alpha): {target_quantile_level:.3f}")

    # The term inside the ceiling function
    term = N * target_quantile_level
    print(f"(n+1) * (1-alpha): {term:.3f}")

    # The rank k, based on the ceiling of (n+1)(1-alpha)
    k = math.ceil(term)
    print(f"Quantile rank k = ceil({term:.3f}): {k}")

    # The lowest possible coverage probability is k / (n+1)
    lowest_coverage = k / N
    print(f"\nFinal Equation: {k} / {N}")
    
    print(f"\nLowest guaranteed coverage probability: {lowest_coverage:.4f}")
    return lowest_coverage

if __name__ == '__main__':
    # --- User-configurable values ---
    # n is the number of training samples
    n_samples = 20
    # alpha is the significance level
    alpha_level = 0.1
    # ----------------------------------
    
    calculate_lowest_coverage(n_samples, alpha_level)
