import math

def calculate_lowest_coverage(n, alpha):
    """
    Calculates the lowest possible coverage probability for a conformal prediction
    interval based on leave-one-out (LOO) residuals. This guarantee holds for any
    data distribution and any algorithm A.

    Args:
        n (int): The number of training data points.
        alpha (float): The desired significance level (e.g., 0.1 for a nominal 90% interval).

    Returns:
        float: The minimum guaranteed coverage probability.
    """
    if not isinstance(n, int) or n <= 0:
        raise ValueError("n must be a positive integer.")
    if not isinstance(alpha, float) or not (0 < alpha < 1):
        raise ValueError("alpha must be a float between 0 and 1.")

    # The total number of exchangeable residuals (n training points + 1 test point).
    m = n + 1

    # The quantile level k is derived from the target coverage 1-alpha.
    # We use ceiling to ensure the property holds.
    # k = ceil(m * (1-alpha))
    k_numerator = math.ceil(m * (1 - alpha))

    # The lowest coverage probability is guaranteed to be k / m.
    min_coverage = k_numerator / m
    
    # Per the instructions, we output each number in the final equation.
    print(f"For a sample size of n = {n} and a significance level of alpha = {alpha}:")
    print("The lowest guaranteed coverage probability is given by the formula:")
    print("  P(coverage) = ceil((n + 1) * (1 - alpha)) / (n + 1)")
    print("")
    print("Plugging in the values:")
    print(f"  P(coverage) = ceil(({n} + 1) * (1 - {alpha})) / ({n} + 1)")
    print(f"              = ceil({m} * {1 - alpha}) / {m}")
    print(f"              = ceil({m * (1-alpha):.2f}) / {m}")
    print(f"              = {k_numerator} / {m}")
    print(f"              = {min_coverage}")
    
    return min_coverage

# --- Example Calculations ---
print("--- Example 1: A standard setting ---")
# Let's consider a case with n=99 training points and alpha=0.1
calculate_lowest_coverage(n=99, alpha=0.1)
# Here, (99+1)*(1-0.1) = 100*0.9 = 90. So ceil(90)=90. The coverage is 90/100 = 0.9.

print("\n" + "="*40 + "\n")

print("--- Example 2: Illustrating the ceiling function ---")
# Let's consider a case where (n+1)*alpha is not an integer.
# Let n=98, alpha=0.1.
# Then (98+1)*(1-0.1) = 99*0.9 = 89.1. ceil(89.1)=90. The coverage is 90/99.
calculate_lowest_coverage(n=98, alpha=0.1)
