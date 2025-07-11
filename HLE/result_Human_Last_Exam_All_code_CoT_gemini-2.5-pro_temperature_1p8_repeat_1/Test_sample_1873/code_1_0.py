import math

def calculate_sampling_ratio(completeness_tolerance, confidence_level, power_law_exponent, pareto_shape):
    """
    Calculates the minimum sampling ratio r based on a model incorporating
    statistical requirements and graph structural properties.

    The formula used is:
    r = ( (alpha - 1) / (gamma - 1) ) * ( epsilon / ln(1 / delta) )
    where delta is 1 - confidence_level.
    """
    
    # Parameters from the problem statement
    epsilon = completeness_tolerance
    gamma = power_law_exponent
    alpha = pareto_shape
    delta = 1 - confidence_level

    # Structural complexity factor from graph properties
    # This represents the ratio of parameters related to the means of the distributions
    structural_complexity = (alpha - 1) / (gamma - 1)

    # Statistical effort term based on tolerance and confidence
    # ln(1/delta) is used based on common bounds in learning theory
    statistical_effort = epsilon / math.log(1 / delta)

    # Calculate the final minimum sampling ratio r
    r = structural_complexity * statistical_effort

    # Print the values used in the calculation
    print(f"Given parameters:")
    print(f"Marginal completeness tolerance (ε): {epsilon}")
    print(f"Confidence level: {confidence_level} (δ = {delta})")
    print(f"Power-law exponent (γ): {gamma}")
    print(f"Pareto distribution shape (α): {alpha}")
    print("-" * 30)
    print("Calculation steps:")
    print(f"Structural Complexity Term ((α - 1) / (γ - 1)) = (({alpha} - 1) / ({gamma} - 1)) = {structural_complexity:.4f}")
    print(f"Statistical Effort Term (ε / ln(1/δ)) = ({epsilon} / ln(1/{delta})) = {statistical_effort:.4f}")
    print(f"Final Ratio (r) = Structural Complexity * Statistical Effort = {structural_complexity:.4f} * {statistical_effort:.4f}")
    print("-" * 30)
    
    # Round the final result to 4 decimal places as requested
    r_rounded = round(r, 4)
    
    print(f"The minimum required sampling ratio (r) is: {r_rounded}")
    return r_rounded

if __name__ == "__main__":
    # Given parameters
    epsilon = 0.05
    confidence = 0.99
    gamma = 2.1
    alpha = 2.5

    # Calculate and get the final answer
    final_ratio = calculate_sampling_ratio(epsilon, confidence, gamma, alpha)
