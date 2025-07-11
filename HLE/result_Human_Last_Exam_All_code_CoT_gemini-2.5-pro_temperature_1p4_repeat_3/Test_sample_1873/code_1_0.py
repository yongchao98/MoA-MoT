import scipy.stats

def calculate_sampling_ratio():
    """
    Calculates the minimum sampling ratio r based on knowledge graph properties.
    """
    # Given parameters
    gamma = 2.1  # Power-law exponent
    alpha = 2.5  # Pareto distribution shape
    epsilon = 0.05  # Marginal completeness tolerance
    confidence_level = 0.99

    # Calculate the Z-score for the given confidence level
    # For a two-tailed interval, we use 1 - (1 - confidence_level) / 2
    z_score = scipy.stats.norm.ppf(1 - (1 - confidence_level) / 2)

    # The formula used to calculate the ratio r
    # r = (Z * ε) / ((γ - 1)² * (α - 1))
    
    # Calculate numerator and denominator
    numerator = z_score * epsilon
    gamma_factor = (gamma - 1)**2
    alpha_factor = (alpha - 1)
    denominator = gamma_factor * alpha_factor
    
    # Calculate the final ratio
    r = numerator / denominator

    # --- Outputting the results ---
    
    # Print the explanation and the formula
    print("The formula to calculate the sampling ratio r is:")
    print("r = (Z * ε) / ((γ - 1)² * (α - 1))")
    print("\nWhere:")
    print(f"  γ (power-law exponent) = {gamma}")
    print(f"  α (Pareto shape) = {alpha}")
    print(f"  ε (marginal tolerance) = {epsilon}")
    print(f"  Z (Z-score for {confidence_level*100}% confidence) = {z_score:.4f}")

    # Print the calculation with substituted values
    print("\nCalculation:")
    print(f"r = ({z_score:.4f} * {epsilon}) / (({gamma} - 1)² * ({alpha} - 1))")
    print(f"r = {numerator:.4f} / (({gamma - 1:.1f})² * ({alpha - 1:.1f}))")
    print(f"r = {numerator:.4f} / ({gamma_factor:.2f} * {alpha_factor:.1f})")
    print(f"r = {numerator:.4f} / {denominator:.3f}")
    
    # Print the final result
    print(f"\nThe calculated minimum sampling ratio r is: {r:.4f}")
    
    # Return the final value for the answer tag
    return round(r, 4)

# Execute the function and capture the final answer
final_answer = calculate_sampling_ratio()