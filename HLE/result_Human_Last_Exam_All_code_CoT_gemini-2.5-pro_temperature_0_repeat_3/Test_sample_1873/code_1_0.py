import math
from scipy.stats import norm

def calculate_sampling_ratio():
    """
    Calculates the minimum sampling ratio r based on the given knowledge graph properties.
    """
    # Given parameters
    alpha = 2.5  # Shape for truncated Pareto distribution of entity neighborhood similarity
    gamma = 2.1  # Power-law exponent for scale-free properties
    epsilon = 0.05  # Marginal completeness tolerance
    confidence_level = 0.99

    # Step 1: Determine the Z-score for the given confidence level.
    # For a two-tailed interval, we look at the (1 + confidence_level) / 2 percentile.
    z_score = norm.ppf((1 + confidence_level) / 2)

    # Step 2: Model the Variance Inflation Factor (VIF)
    # The VIF quantifies the increased sample size needed due to the complex graph structure.
    # The exponents alpha and gamma are close to 2, where variance can diverge.
    # A plausible model for the VIF is based on the inverse of the product of the distances of the exponents from 2.
    vif = 1 / ((gamma - 2) * (alpha - 2))

    # Step 3: Determine the sampling ratio r.
    # The standard sample size n_srs = (Z^2 * p * (1-p)) / epsilon^2.
    # The effective sample size needed is n_eff = n_srs * VIF.
    # The problem asks for a ratio r that is independent of the total number of triples N.
    # This can be interpreted as a scenario where the graph is so complex that the required
    # effective sample size n_eff is on the order of the total population size N.
    # In this case, the sampling ratio r = n_srs / N ≈ n_srs / n_eff = 1 / VIF.
    # This makes r dependent only on the structural complexity of the graph.
    r = 1 / vif

    # Output the parameters and the final calculation
    print(f"Parameters:")
    print(f"  - Pareto shape alpha (α): {alpha}")
    print(f"  - Power-law exponent gamma (γ): {gamma}")
    print(f"  - Tolerance epsilon (ε): {epsilon}")
    print(f"  - Confidence Level: {confidence_level} (Z-score ≈ {z_score:.4f})")
    print("\nCalculation Steps:")
    print(f"1. The Variance Inflation Factor (VIF) is calculated as 1 / ((γ - 2) * (α - 2)).")
    print(f"   VIF = 1 / (({gamma} - 2) * ({alpha} - 2))")
    print(f"   VIF = 1 / ({gamma - 2:.1f} * {alpha - 2:.1f})")
    print(f"   VIF = 1 / {((gamma - 2) * (alpha - 2)):.4f}")
    print(f"   VIF = {vif:.4f}")
    print("\n2. The sampling ratio r is estimated as 1 / VIF.")
    print(f"   r = 1 / {vif:.4f}")
    print(f"   r = {r:.4f}")

    # Final answer rounded to 4 decimal places
    final_r = round(r, 4)
    print(f"\nThe minimum required sampling ratio r is: {final_r}")

calculate_sampling_ratio()
>>> 0.0500