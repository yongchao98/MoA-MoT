import math

def calculate_loo_conformal_coverage(n, alpha):
    """
    Calculates the minimum guaranteed coverage probability for LOO conformal prediction.

    Args:
        n (int): The number of training data points.
        alpha (float): The desired significance level (e.g., 0.1 for 90% confidence).

    Returns:
        float: The minimum coverage probability.
    """
    if not (isinstance(n, int) and n > 0):
        raise ValueError("n must be a positive integer.")
    if not (0 < alpha < 1):
        raise ValueError("alpha must be between 0 and 1.")

    # The total number of items in the exchangeable sequence is n+1
    N = n + 1
    
    # The desired quantile level for the interval
    quantile_level = 1 - alpha
    
    # Calculate k, the rank for the quantile
    # k = ceil((1-alpha) * (n+1))
    k = math.ceil(quantile_level * N)
    
    # The lowest guaranteed coverage probability is k / (n+1)
    min_coverage_prob = k / N
    
    # The problem asks to output the numbers in the equation
    print(f"For n = {n} and alpha = {alpha}:")
    print("The minimum coverage probability is given by the formula: ceil((1 - alpha) * (n + 1)) / (n + 1)")
    print(f"Plugging in the numbers:")
    print(f"n + 1 = {N}")
    print(f"1 - alpha = {quantile_level}")
    print(f"k = ceil({quantile_level} * {N}) = ceil({quantile_level * N}) = {k}")
    print(f"Lowest coverage probability = {k} / {N} = {min_coverage_prob}")
    
    return min_coverage_prob

# Example usage:
n_example = 19
alpha_example = 0.1
calculate_loo_conformal_coverage(n_example, alpha_example)

print("-" * 20)

# Another example:
n_example_2 = 20
alpha_example_2 = 0.05
calculate_loo_conformal_coverage(n_example_2, alpha_example_2)
