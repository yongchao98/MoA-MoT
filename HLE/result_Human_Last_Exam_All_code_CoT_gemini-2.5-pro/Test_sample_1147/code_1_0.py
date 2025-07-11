import math

def calculate_lowest_coverage(n, alpha):
    """
    Calculates the guaranteed lowest coverage probability for a Jackknife+ prediction interval.

    Args:
        n (int): The number of training data points.
        alpha (float): The desired significance level (e.g., 0.1 for 90% target coverage).
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.")
        return
    if not isinstance(alpha, (int, float)) or not (0 < alpha < 1):
        print("Error: alpha must be between 0 and 1.")
        return

    # --- Calculation ---
    n_plus_1 = n + 1
    one_minus_alpha = 1 - alpha
    
    # This is the term (1-alpha)*(n+1)
    product = one_minus_alpha * n_plus_1
    
    # This is the numerator: ceil((1-alpha)*(n+1))
    k = math.ceil(product)
    
    # This is the final guaranteed lowest probability
    lowest_prob = k / n_plus_1

    # --- Output the result with equation steps ---
    print(f"Given n = {n} and alpha = {alpha}:")
    print("The theoretical lowest coverage probability is given by the formula:")
    print("P_min = ceil((1 - alpha) * (n + 1)) / (n + 1)\n")

    print("Substituting the values into the formula:")
    print(f"P_min = ceil((1 - {alpha}) * ({n} + 1)) / ({n} + 1)")
    print(f"P_min = ceil({one_minus_alpha:.2f} * {n_plus_1}) / {n_plus_1}")
    print(f"P_min = ceil({product:.2f}) / {n_plus_1}")
    print(f"P_min = {k} / {n_plus_1}")
    print(f"P_min = {lowest_prob}")
    print(f"\nThis means the coverage is guaranteed to be at least {lowest_prob:.4f}, or {lowest_prob:.2%}.")

# Example usage with some sample values
# You can change these values to see the result for different scenarios
n_train_points = 20
alpha_level = 0.1
calculate_lowest_coverage(n_train_points, alpha_level)
