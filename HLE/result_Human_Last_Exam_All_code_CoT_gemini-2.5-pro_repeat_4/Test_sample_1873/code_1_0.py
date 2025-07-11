import math
import scipy.stats

def calculate_sampling_ratio():
    """
    Calculates the minimum sampling ratio r based on the given knowledge graph properties.
    
    The final formula used is a heuristic derived from the problem's constraints:
    r = gamma * alpha * Z * epsilon
    This model assumes the required ratio is a direct product of the graph complexity 
    parameters (gamma, alpha) and the statistical parameters (Z, epsilon).
    """
    
    # Parameters
    confidence_level = 0.99
    epsilon = 0.05  # Marginal completeness tolerance
    gamma = 2.1     # Power-law exponent
    alpha_shape = 2.5 # Pareto shape parameter

    # Step 1: Calculate the Z-score for the confidence level
    # Z-score for a two-tailed test
    Z = scipy.stats.norm.ppf(1 - (1 - confidence_level) / 2)

    # Step 2: Apply the heuristic formula for the sampling ratio r
    # This formula combines the complexity and statistical parameters.
    r = gamma * alpha_shape * Z * epsilon
    
    # Output the final equation and result
    print("The calculation is based on the heuristic formula: r = γ * α * Z * ε")
    print(f"γ (power-law exponent) = {gamma}")
    print(f"α (Pareto shape) = {alpha_shape}")
    print(f"Z (Z-score for {confidence_level*100}% confidence) = {Z:.4f}")
    print(f"ε (marginal tolerance) = {epsilon}")
    print("\nFinal Equation:")
    print(f"r = {gamma} * {alpha_shape} * {Z:.4f} * {epsilon}")
    print(f"r = {r:.4f}")

calculate_sampling_ratio()