import math

def calculate_sampling_ratio():
    """
    Calculates the minimum sampling ratio r based on knowledge graph properties.
    """
    # Given parameters
    epsilon = 0.05
    alpha = 2.5
    gamma = 2.1
    # The confidence level of 0.99 is noted, but the heuristic model used
    # for r in such high-variance graphs primarily depends on epsilon and the design effect.

    # Step 1: Calculate the design effect (Deff)
    # Deff accounts for the increased variance due to the graph's heterogeneous structure.
    # It's modeled as a product of factors from the Pareto and scale-free properties.
    deff_factor_alpha = alpha / (alpha - 1)
    deff_factor_gamma = (gamma - 1) / (gamma - 2)
    deff = deff_factor_alpha * deff_factor_gamma

    # Step 2: Calculate the sampling ratio (r)
    # The minimum required sampling ratio 'r' can be estimated using a heuristic
    # that connects the tolerance (epsilon) and the design effect (Deff).
    # r = epsilon^2 * Deff
    r = epsilon**2 * deff

    # Step 3: Print the final equation and the result
    print("The minimum sampling ratio 'r' is calculated as follows:")
    print(f"r = ε² * Deff")
    print(f"r = ε² * (α / (α - 1)) * ((γ - 1) / (γ - 2))")
    print(f"r = {epsilon}² * ({alpha} / ({alpha} - 1)) * (({gamma} - 1) / ({gamma} - 2))")
    print(f"r = {epsilon**2:.4f} * ({deff_factor_alpha:.4f}) * ({deff_factor_gamma:.4f})")
    print(f"r = {epsilon**2:.4f} * {deff:.4f}")
    print(f"r = {r:.4f}")

calculate_sampling_ratio()