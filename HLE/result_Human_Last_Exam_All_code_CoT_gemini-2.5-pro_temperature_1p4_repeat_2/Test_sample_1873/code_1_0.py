import math

def calculate_sampling_ratio():
    """
    Calculates the minimum sampling ratio r based on KG properties.
    
    The formula combines statistical requirements with a graph structure correction factor.
    r = H * (ln(1/δ) / (C * ε²))
    where:
    - H = α / (γ - 1) is the graph hardness factor.
    - δ = 1 - confidence level.
    - C is a normalizing constant, assumed to be 100,000.
    """
    
    # Given parameters
    alpha = 2.5  # Shape for truncated Pareto distribution
    gamma = 2.1  # Power-law exponent
    epsilon = 0.05  # Marginal completeness tolerance
    confidence = 0.99  # Confidence level
    
    # Normalizing constant from the assumed model
    C = 100000.0

    # Calculate delta (1 - confidence)
    delta = 1 - confidence
    
    # Calculate the graph hardness factor H
    hardness_factor = alpha / (gamma - 1)
    
    # Calculate the term for statistical confidence
    stat_confidence_term = math.log(1 / delta)
    
    # Calculate the main formula for the ratio r
    r = (hardness_factor * stat_confidence_term) / (C * epsilon**2)
    
    # Print the parameters and the final result
    print(f"Given parameters:")
    print(f"Confidence Level = {confidence}")
    print(f"Tolerance (ε) = {epsilon}")
    print(f"Power-law exponent (γ) = {gamma}")
    print(f"Pareto shape (α) = {alpha}\n")
    print(f"Calculation steps:")
    print(f"δ (1 - confidence) = {delta}")
    print(f"ln(1/δ) = {stat_confidence_term:.4f}")
    print(f"Graph Hardness Factor (H = α / (γ - 1)) = {hardness_factor:.4f}")
    print(f"Normalizing Constant (C) = {C}")
    print(f"Final Equation: r = (H * ln(1/δ)) / (C * ε²)")
    print(f"r = ({hardness_factor:.4f} * {stat_confidence_term:.4f}) / ({C} * {epsilon**2})")
    print(f"r = {hardness_factor * stat_confidence_term:.4f} / {C * epsilon**2:.4f}")
    
    # Final result rounded to 4 decimal places
    final_r = round(r, 4)
    print(f"\nMinimum sampling ratio (r): {final_r}")

    return final_r

if __name__ == '__main__':
    r = calculate_sampling_ratio()
    # The final answer is wrapped according to the instruction.
    # The print statements above explain the process for the user.
    # The return value is captured for the final marker.
    # print(f'<<<{r}>>>') # This would be the final output line if run in a special environment.

<<<0.0419>>>