import scipy.stats

def calculate_sampling_ratio():
    """
    Calculates the minimum sampling ratio r based on the given knowledge graph properties.
    """
    # Given parameters
    confidence_level = 0.99
    epsilon = 0.05  # marginal completeness tolerance
    gamma = 2.1     # power-law exponent
    alpha = 2.5     # Pareto distribution shape

    # Step 1: Find the Z-score for the confidence level
    # The confidence level is 0.99, so the significance level is 0.01.
    # For a two-tailed interval, we use alpha_err/2.
    alpha_err = 1 - confidence_level
    Z = scipy.stats.norm.ppf(1 - alpha_err / 2)

    # Step 2: Calculate the Design Effect (Deff)
    # This captures the extra sample size needed due to the graph's complex structure.
    # Deff = α / (γ - 2)
    deff = alpha / (gamma - 2)

    # Step 3: Calculate the statistical precision factor
    # This is the squared ratio of Z to epsilon.
    stat_precision_factor_sq = (Z / epsilon)**2

    # Step 4: Calculate the final sampling ratio r
    # r = Deff / (Z/ε)²
    r = deff / stat_precision_factor_sq

    # Print the values used in the calculation
    print(f"Confidence Level: {confidence_level}")
    print(f"Z-score (Z): {Z:.4f}")
    print(f"Marginal Tolerance (ε): {epsilon}")
    print(f"Power-law Exponent (γ): {gamma}")
    print(f"Pareto Shape (α): {alpha}")
    print("-" * 30)
    print(f"Design Effect (Deff = α / (γ - 2)): {deff:.4f}")
    print(f"Statistical Precision Factor ((Z/ε)²): {stat_precision_factor_sq:.4f}")
    print("-" * 30)
    print(f"Final Equation: r = Deff / (Z/ε)²")
    print(f"r = {deff:.4f} / ({Z:.4f} / {epsilon})²")
    print(f"r = {deff:.4f} / {stat_precision_factor_sq:.4f}")
    print(f"Calculated sampling ratio r: {r:.4f}")

if __name__ == "__main__":
    calculate_sampling_ratio()
