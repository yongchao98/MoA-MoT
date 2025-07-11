import math
from scipy.stats import norm

def solve():
    """
    Calculates the minimum sampling ratio r based on the given knowledge graph properties.
    """
    # 1. Define the given parameters
    confidence_level = 0.99
    epsilon = 0.05
    p = 0.5  # Assumed for worst-case variance in proportion estimation
    gamma = 2.1 # Power-law exponent for scale-free property
    alpha = 2.5 # Shape for truncated Pareto distribution of entity neighborhood similarity

    # 2. Calculate the Z-score for the given confidence level
    # delta is the significance level (1 - confidence_level)
    delta = 1 - confidence_level
    # Z-score for a two-tailed confidence interval
    z_score = norm.ppf(1 - delta / 2)

    # 3. Calculate the sample size requirement for a simple random sample (n_srs)
    # This represents the "sampling difficulty"
    n_srs = (z_score**2 * p * (1 - p)) / (epsilon**2)
    
    # 4. Calculate the Design Effect (Deff) from graph heterogeneity parameters
    # This represents the "complexity difficulty" due to network structure
    # This factor models the inflation in variance from heavy-tailed distributions
    deff_gamma = (gamma - 1) / (gamma - 2)
    deff_alpha = (alpha - 1) / (alpha - 2)
    deff = deff_gamma * deff_alpha

    # 5. Calculate the ratio r using a heuristic formula
    # Standard formula r = n/N is not usable as N is unknown.
    # This formula combines the sampling difficulty and complexity difficulty into a ratio.
    r = deff / (n_srs + deff)

    # Print the values used in the calculation
    print(f"Confidence Level: {confidence_level}")
    print(f"Z-score (Z): {z_score:.4f}")
    print(f"Marginal Tolerance (ε): {epsilon}")
    print(f"Power-law exponent (γ): {gamma}")
    print(f"Pareto shape (α): {alpha}")
    print("----")
    print(f"Calculated base sample size (n_srs): {n_srs:.4f}")
    print(f"Calculated design effect (Deff): {deff:.4f}")
    print("----")
    # Final equation showing how r is computed
    print(f"Final ratio r = Deff / (n_srs + Deff)")
    print(f"r = {deff:.4f} / ({n_srs:.4f} + {deff:.4f})")
    print(f"r = {deff:.4f} / {n_srs + deff:.4f}")
    print(f"r = {r:.4f}")

    # Final result rounded to 4 decimal places
    final_r = round(r, 4)
    print(f"\nThe minimum required sampling ratio r is: {final_r}")

solve()
<<<0.0474>>>