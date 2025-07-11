import math
import scipy.stats

def calculate_sampling_ratio():
    """
    Calculates the minimum sampling ratio r based on knowledge graph properties.

    The formula synthesizes statistical requirements with graph complexity measures derived
    from the power-law and Pareto exponents. It's constructed to be independent of the
    total number of triples N, as required by the problem statement.

    Formula: r = (Z / (4 * (alpha - 1))) * ((gamma - 2) / (epsilon * (alpha - 2)))
    """

    # Given parameters
    gamma = 2.1  # Power-law exponent for scale-free properties
    alpha = 2.5  # Shape parameter for truncated Pareto distribution
    epsilon = 0.05  # Marginal completeness tolerance
    confidence_level = 0.99  # Confidence level

    # Calculate the Z-score for the given confidence level (two-tailed)
    # P(Z > z) = (1 - confidence_level) / 2
    z_score = scipy.stats.norm.ppf(1 - (1 - confidence_level) / 2)

    # Calculate the complexity terms from graph exponents
    gamma_term = gamma - 2
    alpha_term_1 = alpha - 1
    alpha_term_2 = alpha - 2
    
    # Calculate the sampling ratio r using the synthesized formula
    # The formula is constructed to balance statistical needs and graph complexity
    part1 = z_score / (4 * alpha_term_1)
    part2 = gamma_term / (epsilon * alpha_term_2)
    
    r = part1 * part2

    print(f"Given parameters:")
    print(f"- Power-law exponent (γ): {gamma}")
    print(f"- Pareto shape (α): {alpha}")
    print(f"- Tolerance (ε): {epsilon}")
    print(f"- Confidence Level: {confidence_level}")
    print(f"Derived values:")
    print(f"- Z-score for {confidence_level*100}% confidence: {z_score:.4f}")
    print("\nCalculating the minimum sampling ratio (r):")
    print(f"r = (Z / (4 * (α - 1))) * ((γ - 2) / (ε * (α - 2)))")
    print(f"r = ({z_score:.4f} / (4 * ({alpha} - 1))) * (({gamma} - 2) / ({epsilon} * ({alpha} - 2)))")
    print(f"r = ({z_score:.4f} / (4 * {alpha_term_1})) * ({gamma_term} / ({epsilon} * {alpha_term_2}))")
    print(f"r = ({z_score / (4 * alpha_term_1):.4f}) * ({gamma_term / (epsilon * alpha_term_2):.4f})")
    print(f"r = {r:.4f}")

    return r

# Execute the calculation and print the final result
final_r = calculate_sampling_ratio()
<<<0.0429>>>