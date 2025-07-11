import math

def calculate_loo_coverage(n, alpha):
    """
    Calculates the exact theoretical coverage probability for Leave-One-Out
    conformal prediction.

    Args:
        n (int): The number of training samples.
        alpha (float): The desired significance level (e.g., 0.1 for 90% confidence).
    
    Returns:
        float: The guaranteed coverage probability.
    """
    if not (isinstance(n, int) and n > 0):
        raise ValueError("n must be a positive integer.")
    if not (isinstance(alpha, float) and 0 < alpha < 1):
        raise ValueError("alpha must be a float between 0 and 1.")

    n_plus_1 = n + 1
    one_minus_alpha = 1 - alpha
    
    # Term inside the ceiling function
    term = n_plus_1 * one_minus_alpha
    
    # The numerator of the coverage formula is k
    k = math.ceil(term)
    
    # The exact coverage probability
    coverage = k / n_plus_1
    
    # Print the step-by-step calculation
    print(f"For n = {n} and alpha = {alpha}:")
    print("The coverage probability P is given by ceil((n+1)*(1-alpha)) / (n+1)")
    print(f"P = ceil(({n}+1) * (1-{alpha})) / ({n}+1)")
    print(f"P = ceil({n_plus_1} * {one_minus_alpha}) / {n_plus_1}")
    print(f"P = ceil({term:.4f}) / {n_plus_1}")
    print(f"P = {k} / {n_plus_1}")
    print(f"Guaranteed Coverage = {coverage:.4f}")
    
    # Compare with the 1 - alpha bound
    print(f"\nThis value is always >= (1 - alpha), which is {one_minus_alpha:.4f}.")
    print("-" * 20)

# Example 1: (n+1)(1-alpha) is not an integer
calculate_loo_coverage(n=10, alpha=0.1)

# Example 2: (n+1)(1-alpha) is an integer
calculate_loo_coverage(n=99, alpha=0.1)

# Example 3: alpha is small enough that coverage becomes 1
calculate_loo_coverage(n=20, alpha=0.04)
