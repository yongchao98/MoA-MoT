import scipy.stats
import math

def calculate_sampling_ratio():
    """
    Calculates the minimum sampling ratio r based on the given knowledge graph properties.
    """
    # Step 1: Define the given parameters
    alpha = 2.5  # Shape parameter for Pareto distribution
    gamma = 2.1  # Power-law exponent for scale-free properties
    epsilon = 0.05  # Marginal completeness tolerance
    confidence_level = 0.99

    # Step 2: Calculate the Z-score for the given confidence level
    # The Z-score corresponds to the (1 - (1 - confidence_level)/2) percentile
    z_score = scipy.stats.norm.ppf(1 - (1 - confidence_level) / 2)

    # Step 3: Assume worst-case completeness (p=0.5) for maximum variance
    p = 0.5
    p_variance = p * (1 - p)

    # Step 4: Calculate the required sample size for a Simple Random Sample (SRS)
    # This represents the statistical sample requirement.
    n0 = (z_score**2 * p_variance) / (epsilon**2)

    # Step 5: Calculate the Design Effect (Deff) from the graph's structural properties
    # The heterogeneity factor for a power-law-like distribution with exponent x is (x-1)/(x-2).
    # This factor is related to the ratio of the first moment to the minimum value.
    deff_gamma = (gamma - 1) / (gamma - 2)
    deff_alpha = (alpha - 1) / (alpha - 2)

    # The total design effect is the product of the individual effects.
    total_deff = deff_gamma * deff_alpha

    # Step 6: Calculate the final ratio r.
    # Based on the interpretation of the problem, the ratio r is derived from the
    # relationship between the structural complexity (Deff) and the statistical
    # sample requirement (n0). A plausible model is r = Deff / n0.
    r = total_deff / n0
    
    # Output the intermediate values to show the calculation steps
    print(f"Confidence Level: {confidence_level}")
    print(f"Z-score: {z_score:.4f}")
    print(f"Margin of Error (epsilon): {epsilon}")
    print(f"Power-law exponent (gamma): {gamma}")
    print(f"Pareto shape (alpha): {alpha}")
    print(f"---")
    print(f"Statistical Sample Requirement (n0): (Z^2 * p(1-p)) / epsilon^2 = ({z_score:.4f}^2 * {p_variance}) / {epsilon}^2 = {n0:.4f}")
    print(f"Design Effect from gamma: (gamma-1)/(gamma-2) = ({gamma-1})/({gamma-2}) = {deff_gamma:.4f}")
    print(f"Design Effect from alpha: (alpha-1)/(alpha-2) = ({alpha-1})/({alpha-2}) = {deff_alpha:.4f}")
    print(f"Total Design Effect (Deff): {deff_gamma:.4f} * {deff_alpha:.4f} = {total_deff:.4f}")
    print(f"---")
    print(f"Final Ratio (r) = Deff / n0 = {total_deff:.4f} / {n0:.4f} = {r:.4f}")
    print(f"\nThe minimum required sampling ratio r is {r:.4f}.")

calculate_sampling_ratio()