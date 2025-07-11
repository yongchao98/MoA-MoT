import scipy.stats as st

def calculate_sampling_ratio():
    """
    Calculates the minimum sampling ratio r based on the given knowledge graph properties.
    """
    # Parameters from the problem statement
    confidence_level = 0.99
    epsilon = 0.05  # Marginal completeness tolerance
    gamma = 2.1     # Scale-free power-law exponent
    alpha = 2.5     # Pareto distribution shape

    # Step 1: Calculate the Z-score for the given confidence level
    # For a two-tailed test, alpha_err is (1 - confidence_level)
    # We look up the Z-score for 1 - alpha_err / 2
    z_score = st.norm.ppf(1 - (1 - confidence_level) / 2)

    # Step 2: Apply the formula for the sampling ratio r
    # r = (Z * ε) / ((γ - 1) * (α - 1))
    numerator = z_score * epsilon
    denominator = (gamma - 1) * (alpha - 1)
    r = numerator / denominator

    # Step 3: Round the result to 4 decimal places
    r_rounded = round(r, 4)

    # Print the explanation and the final equation
    print("The minimum sampling ratio r is calculated using the formula:")
    print("r = (Z * ε) / ((γ - 1) * (α - 1))")
    print("\nGiven values:")
    print(f"Confidence Level = {confidence_level} -> Z-score = {z_score:.4f}")
    print(f"Marginal Tolerance (ε) = {epsilon}")
    print(f"Scale-free Exponent (γ) = {gamma}")
    print(f"Pareto Shape (α) = {alpha}")
    print("\nCalculation:")
    print(f"r = ({z_score:.4f} * {epsilon}) / (({gamma} - 1) * ({alpha} - 1))")
    print(f"r = {numerator:.4f} / (({gamma - 1:.1f}) * ({alpha - 1:.1f}))")
    print(f"r = {numerator:.4f} / {denominator:.2f}")
    print(f"r ≈ {r:.7f}")
    print(f"\nRounding to 4 decimal places, the minimum sampling ratio r is: {r_rounded}")

calculate_sampling_ratio()
<<<0.0781>>>