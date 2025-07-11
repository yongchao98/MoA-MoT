import math

def calculate_sampling_ratio():
    """
    Calculates the minimum sampling ratio 'r' required to estimate predicate
    completeness in a knowledge graph with given properties.
    """
    # Given parameters from the problem description
    gamma = 2.1  # Power-law exponent of the scale-free graph
    alpha = 2.5  # Shape of the truncated Pareto distribution for neighborhood similarity
    epsilon = 0.05  # Marginal completeness tolerance
    confidence = 0.99  # Desired confidence level

    # Calculate delta (1 - confidence)
    delta = 1 - confidence

    # The formula used to calculate the sampling ratio 'r' is:
    # r = ln(1/δ) / ( (1/(γ-2)) * (1/ε) * (α-1) )
    
    # Calculate each component of the formula
    log_term = math.log(1 / delta)
    gamma_factor = 1 / (gamma - 2)
    epsilon_factor = 1 / epsilon
    alpha_factor = alpha - 1
    
    # Calculate the full denominator
    denominator = gamma_factor * epsilon_factor * alpha_factor
    
    # Calculate the final ratio 'r'
    r = log_term / denominator
    
    # Print the explanation and step-by-step calculation
    print("To find the minimum required sampling ratio 'r', we use the following formula:")
    print("r = ln(1/δ) / ( (1/(γ-2)) * (1/ε) * (α-1) )")
    print("\n--- Given Values ---")
    print(f"Power-law exponent (γ) = {gamma}")
    print(f"Pareto shape (α) = {alpha}")
    print(f"Tolerance (ε) = {epsilon}")
    print(f"Confidence = {confidence} which means δ = {delta:.2f}")
    
    print("\n--- Calculation ---")
    print(f"r = ln(1/{delta:.2f}) / ( (1/({gamma} - 2)) * (1/{epsilon}) * ({alpha} - 1) )")
    print(f"r = {log_term:.4f} / ( (1/{gamma - 2:.1f}) * ({epsilon_factor}) * ({alpha_factor:.1f}) )")
    print(f"r = {log_term:.4f} / ( {gamma_factor:.1f} * {epsilon_factor} * {alpha_factor:.1f} )")
    print(f"r = {log_term:.4f} / {denominator:.1f}")
    print(f"The calculated ratio is r = {r}")
    print(f"\nFinal Answer (rounded to 4 decimal places): {r:.4f}")

# Execute the function to get the result
calculate_sampling_ratio()
<<<0.0154>>>