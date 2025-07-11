import math
from scipy.stats import norm

def solve():
    """
    Calculates the minimum sampling ratio r required to estimate predicate completeness.
    """
    # Given parameters
    alpha = 2.5
    gamma = 2.1
    epsilon = 0.05
    confidence_level = 0.99

    # --- Step 1: Calculate the Z-score for the confidence level ---
    # The Z-score for a two-tailed confidence interval.
    z_score = norm.ppf(1 - (1 - confidence_level) / 2)

    # --- Step 2: Determine worst-case proportion for maximum variance ---
    # The variance of a proportion is maximized when p = 0.5.
    p = 0.5

    # --- Step 3: Calculate the base sample size (n_base) ---
    # This is the standard sample size for a simple random sample.
    n_base = (z_score**2 * p * (1 - p)) / (epsilon**2)

    # --- Step 4: Calculate the Structural Complexity Factor (C) ---
    # This factor models the variance inflation due to the graph's heavy-tailed distributions.
    # The term (3-exponent) is characteristic for phenomena in the 2 < exponent < 3 regime.
    try:
        complexity_factor = 1 / ((3 - alpha) * (3 - gamma))
    except ZeroDivisionError:
        print("Error: alpha or gamma cannot be 3 for this formula.")
        return

    # --- Step 5: Calculate the final effective sample size ---
    # The actual required sample size is the base size multiplied by the complexity factor.
    n_final = n_base * complexity_factor
    
    # --- Step 6: Calculate the ratio r ---
    # The ratio 'r' is interpreted as the inverse of the required sample count (n_final).
    # This provides a dimensionless value that depends on all parameters in a consistent way
    # (higher complexity or required precision leads to a smaller r, indicating
    # that each sample contributes less).
    r = 1 / n_final

    # --- Final Output ---
    # Print the full equation with intermediate values.
    print(f"The minimum required sampling ratio r is calculated as follows:")
    print(f"r = 1 / ( (Z^2 * p * (1-p) / ε^2) * (1 / ((3-α)*(3-γ))) )")
    print(f"r = 1 / ( ({z_score:.4f}^2 * {p} * (1-{p}) / {epsilon}^2) * (1 / ((3-{alpha})*(3-{gamma}))) )")
    print(f"r = 1 / ( ({n_base:.4f}) * ({complexity_factor:.4f}) )")
    print(f"r = 1 / {n_final:.4f}")
    print(f"r = {r:.4f}")

solve()
<<<0.0007>>>