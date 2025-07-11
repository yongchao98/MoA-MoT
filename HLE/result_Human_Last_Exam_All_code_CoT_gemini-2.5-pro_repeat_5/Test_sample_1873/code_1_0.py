import scipy.stats

def calculate_sampling_ratio():
    """
    Calculates the minimum sampling ratio r based on the given knowledge graph properties.
    """
    # 1. Define the given parameters
    gamma = 2.1  # Power-law exponent
    alpha = 2.5  # Pareto shape for entity neighborhood similarity
    epsilon = 0.05  # Marginal completeness tolerance
    confidence_level = 0.99

    # 2. Calculate the Z-score for the given confidence level
    # For a two-tailed interval, we use (1 - confidence_level) / 2
    alpha_confidence = 1 - confidence_level
    z_score = scipy.stats.norm.ppf(1 - alpha_confidence / 2)

    # 3. Apply the formula to calculate the sampling ratio r
    # The model assumes r = structural_factor / statistical_requirement
    structural_factor = (gamma - 2) * (alpha - 1)
    statistical_requirement = z_score / epsilon
    
    r = structural_factor / statistical_requirement
    
    # 4. Print the final equation and the result
    print("The calculation is based on the formula: r = ((γ - 2) * (α - 1)) / (Z / ε)")
    print("\nPlugging in the values:")
    print(f"r = (({gamma} - 2) * ({alpha} - 1)) / ({z_score:.4f} / {epsilon})")
    print(f"r = (({gamma - 2:.1f}) * ({alpha - 1:.1f})) / ({statistical_requirement:.4f})")
    print(f"r = {structural_factor:.2f} / {statistical_requirement:.4f}")

    print(f"\nThe minimum required sampling ratio r is: {r:.4f}")

if __name__ == "__main__":
    calculate_sampling_ratio()