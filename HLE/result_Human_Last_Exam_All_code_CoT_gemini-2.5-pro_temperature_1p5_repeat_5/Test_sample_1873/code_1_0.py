import math

def calculate_sampling_ratio():
    """
    Calculates the minimum sampling ratio 'r' based on the problem description.

    The method uses Chebyshev's inequality for a conservative sample size estimate,
    which is appropriate given the scale-free nature of the graph (γ < 3).
    """

    # Parameters from the user query
    confidence_level = 0.99
    marginal_tolerance_epsilon = 0.05
    power_law_exponent_gamma = 2.1  # Justification for using a conservative bound
    
    # We use the worst-case proportion (p=0.5) to maximize the variance,
    # ensuring the sample size is sufficient.
    proportion_p = 0.5
    
    # The significance level (delta) is 1 minus the confidence level
    delta = 1 - confidence_level
    
    # --- Step 1: Calculate the required sample size 'n' using Chebyshev's inequality ---
    # Formula: n >= p*(1-p) / (delta * epsilon^2)
    numerator = proportion_p * (1 - proportion_p)
    denominator = delta * (marginal_tolerance_epsilon ** 2)
    
    sample_size_n = numerator / denominator

    # --- Step 2: Calculate the ratio 'r' ---
    # The prompt asks for a ratio 'r' but does not provide a total number of triples.
    # A plausible interpretation in this context is r = 1/n.
    if sample_size_n > 0:
        ratio_r = 1 / sample_size_n
    else:
        ratio_r = 0

    # Round the final result to 4 decimal places
    rounded_r = round(ratio_r, 4)

    # Output the logic and the final result
    print("Plan:")
    print("1. Use Chebyshev's inequality for a conservative sample size due to the graph's scale-free nature (γ < 3).")
    print("2. The sample size 'n' is calculated using the formula: n = (p*(1-p)) / (δ * ε²).")
    print("3. Use worst-case proportion p=0.5, with δ = 1-0.99=0.01 and ε=0.05.")
    print("\nCalculation:")
    print(f"p = {proportion_p}")
    print(f"δ = {delta}")
    print(f"ε = {marginal_tolerance_epsilon}")
    print(f"n = ({proportion_p} * (1-{proportion_p})) / ({delta} * {marginal_tolerance_epsilon}²) = {sample_size_n}")
    print(f"4. The ratio 'r' is interpreted as 1/n, since the total number of triples is unknown.")
    print(f"r = 1 / {sample_size_n} = {ratio_r}")
    print(f"\nFinal Answer (rounded to 4 decimal places): {rounded_r}")

calculate_sampling_ratio()
<<<0.0001>>>