import math
from scipy.stats import norm

def calculate_sampling_ratio():
    """
    Calculates the minimum sampling ratio r based on the provided knowledge graph properties.
    """
    # 1. Given parameters
    gamma = 2.1  # Power-law exponent
    alpha = 2.5  # Pareto shape parameter
    epsilon = 0.05  # Marginal completeness tolerance
    confidence = 0.99  # Confidence level

    # 2. Determine statistical values
    # Significance level (delta)
    delta = 1 - confidence
    # Z-score for a two-tailed test
    Z = norm.ppf(1 - delta / 2)

    # 3. Formulate and calculate the value
    # Based on models for sample size estimation in heterogeneous KGs (e.g., KG-Sizer),
    # a required sample size S_N (independent of stratum size) can be calculated.
    # S_N = Z**2 / (epsilon**2 * (gamma - 1) * (alpha - 1))
    # The problem asks for a ratio 'r'. As the relationship between S_N and r is
    # not specified, and r=S_N/N would depend on N, we choose to compute r = 1/S_N
    # as a plausible interpretation that yields a constant ratio.
    
    # Numerator of the ratio r
    numerator = epsilon**2 * (gamma - 1) * (alpha - 1)
    
    # Denominator of the ratio r
    denominator = Z**2

    # Calculate the ratio r
    r = numerator / denominator

    # 4. Print the final equation and result
    print(f"The minimum sampling ratio r is calculated using the formula:")
    print(f"r = (ε² * (γ - 1) * (α - 1)) / Z²")
    print(f"Given values:")
    print(f"  γ = {gamma}")
    print(f"  α = {alpha}")
    print(f"  ε = {epsilon}")
    print(f"  Confidence = {confidence} -> Z = {Z:.4f}")
    print(f"\nPlugging in the values:")
    print(f"r = ({epsilon}² * ({gamma} - 1) * ({alpha} - 1)) / {Z:.4f}²")
    print(f"r = ({epsilon**2:.4f} * {gamma - 1:.1f} * {alpha - 1:.1f}) / {Z**2:.4f}")
    print(f"r = {numerator:.6f} / {denominator:.4f}")
    print(f"r = {r:.10f}")
    
    # Round the final result to 4 decimal places for the final answer
    r_rounded = round(r, 4)
    print(f"\nThe calculated ratio r rounded to 4 places is: {r_rounded}")
    print(f"\n<<<r = {r:.4f}>>>")

calculate_sampling_ratio()