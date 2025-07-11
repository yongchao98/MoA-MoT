import math

def calculate_coverage(n, alpha):
    """
    Calculates the exact theoretical coverage probability for Jackknife+ (LOO)
    conformal prediction.

    Args:
        n (int): The number of training points.
        alpha (float): The significance level (e.g., 0.1 for 90% confidence).
    """
    if not (isinstance(n, int) and n > 0):
        print("Error: n must be a positive integer.")
        return
    if not (0 < alpha < 1):
        print("Error: alpha must be between 0 and 1.")
        return

    N = n + 1
    one_minus_alpha = 1 - alpha

    # Calculate k, the rank for the quantile
    # k = ceil((n+1)*(1-alpha))
    k = math.ceil(N * one_minus_alpha)

    # The exact coverage probability is k / (n+1)
    exact_coverage = k / N
    
    # The theoretical lower bound is 1 - alpha
    lower_bound = one_minus_alpha

    print(f"--- Calculating for n={n}, alpha={alpha} ---")
    print("The final equation for coverage probability is: ceil((n+1)*(1-alpha)) / (n+1)")
    print(f"Plugging in the numbers: ceil(({n}+1)*(1-{alpha})) / ({n}+1)")
    print(f"= ceil({N*one_minus_alpha:.2f}) / {N}")
    print(f"= {k} / {N}")
    print(f"Exact Coverage Probability: {exact_coverage:.4f}")
    print(f"Theoretical Lower Bound (1-alpha): {lower_bound:.4f}")
    print(f"Is coverage >= lower bound? {'Yes' if exact_coverage >= lower_bound else 'No'}\n")

# Example 1: A case where coverage is strictly greater than 1-alpha
# (n+1)(1-alpha) = 16 * 0.9 = 14.4, which is not an integer.
calculate_coverage(n=15, alpha=0.1)

# Example 2: A case where coverage is exactly 1-alpha
# (n+1)(1-alpha) = 20 * 0.9 = 18, which is an integer.
calculate_coverage(n=19, alpha=0.1)

# Example 3: A different alpha value
# (n+1)(1-alpha) = 51 * 0.95 = 48.45, which is not an integer.
calculate_coverage(n=50, alpha=0.05)