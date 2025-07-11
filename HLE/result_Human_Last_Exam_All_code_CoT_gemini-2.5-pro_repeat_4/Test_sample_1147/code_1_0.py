import math

def calculate_loo_conformal_coverage(n, alpha):
    """
    Calculates the exact theoretical coverage probability for LOO conformal prediction.

    Args:
        n (int): The number of training points.
        alpha (float): The desired significance level (e.g., 0.1 for 90% coverage).
    """
    print(f"--- Calculating for n={n}, alpha={alpha} ---")

    # The number of scores in the calibration set (including +infinity) is n+1
    num_scores = n + 1
    
    # The target quantile level
    quantile_level = 1 - alpha
    
    # The argument inside the ceiling function
    ceil_arg = num_scores * quantile_level
    
    # The numerator of the coverage formula is the ceiling of the argument
    numerator = math.ceil(ceil_arg)
    
    # The denominator is the number of scores
    denominator = num_scores
    
    # The exact coverage probability
    coverage = numerator / denominator
    
    # The lower bound
    lower_bound = 1 - alpha
    
    print("Equation: P(coverage) = ceil((n+1) * (1-alpha)) / (n+1)")
    print(f"Calculation: ceil(({n}+1) * (1-{alpha})) / ({n}+1)")
    print(f"             = ceil({num_scores} * {quantile_level}) / {num_scores}")
    print(f"             = ceil({ceil_arg:.4f}) / {denominator}")
    print(f"             = {numerator} / {denominator}")
    print(f"             = {coverage:.4f}")
    
    print(f"\nThis coverage of {coverage:.4f} is compared to the nominal level 1-alpha = {lower_bound:.4f}.")
    if abs(coverage - lower_bound) < 1e-9:
        print("The coverage is exactly equal to the lower bound 1-alpha.")
    else:
        print("The coverage is strictly greater than the lower bound 1-alpha.")
    print("-" * 35 + "\n")

if __name__ == '__main__':
    # Case 1: (n+1)*(1-alpha) is an integer. The coverage hits the lower bound.
    n1, alpha1 = 19, 0.1
    calculate_loo_conformal_coverage(n1, alpha1)
    
    # Case 2: (n+1)*(1-alpha) is not an integer. The coverage is above the lower bound.
    n2, alpha2 = 18, 0.1
    calculate_loo_conformal_coverage(n2, alpha2)