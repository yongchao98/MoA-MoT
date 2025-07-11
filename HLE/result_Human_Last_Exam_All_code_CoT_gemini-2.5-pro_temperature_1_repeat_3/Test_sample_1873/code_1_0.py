import math

def calculate_sampling_ratio():
    """
    Calculates the minimum sampling ratio r based on the given knowledge graph properties.
    """
    
    # Given parameters
    gamma = 2.1  # Power-law exponent
    alpha = 2.5  # Pareto distribution shape for similarity
    epsilon = 0.05  # Marginal completeness tolerance
    confidence = 0.99  # Confidence level

    # The probability of error
    delta = 1 - confidence

    # Structural complexity factor based on the graph parameters.
    # We use (gamma - 1) and (alpha - 1) as they relate to the first moments
    # of the respective power-law distributions. The inverse represents difficulty.
    structural_difficulty = 1 / ((gamma - 1) * (alpha - 1))

    # Statistical requirement factor.
    # In this model, we assume a direct relationship based on tolerance and error probability.
    statistical_difficulty = delta / epsilon

    # The formula for the sampling ratio r combines these factors.
    r = statistical_difficulty * structural_difficulty
    
    # The final equation with all numbers plugged in
    print(f"The calculation is based on the formula: r = ( (1 - confidence) / epsilon ) * ( 1 / ((gamma - 1) * (alpha - 1)) )")
    print(f"Plugging in the values:")
    print(f"r = ( (1 - {confidence}) / {epsilon} ) * ( 1 / (({gamma} - 1) * ({alpha} - 1)) )")
    print(f"r = ( {delta} / {epsilon} ) * ( 1 / ({gamma-1} * {alpha-1}) )")
    print(f"r = {delta/epsilon} * ( 1 / { (gamma-1)*(alpha-1) } )")
    print(f"r = {delta/epsilon} * {1/((gamma-1)*(alpha-1))}")
    print(f"r = {r}")

    # Round the result to 4 decimal places
    r_rounded = round(r, 4)
    print(f"\nThe minimum required sampling ratio r, rounded to 4 decimal places, is: {r_rounded}")

calculate_sampling_ratio()
<<<0.1212>>>