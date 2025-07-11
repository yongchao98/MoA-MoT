import math
from scipy.stats import norm

def calculate_sampling_ratio():
    """
    Calculates the minimum sampling ratio r based on knowledge graph properties.
    """
    # Given parameters
    gamma = 2.1  # Power-law exponent for scale-free properties
    alpha = 2.5  # Shape for Pareto distribution of neighborhood similarity
    epsilon = 0.05  # Marginal completeness tolerance
    confidence = 0.99 # Confidence level

    # Step 1: Calculate the Z-score for the given confidence level
    # For a two-tailed interval, alpha_risk is (1 - confidence)
    # We look up the Z-score for confidence + alpha_risk / 2
    z_score = norm.ppf(1 - (1 - confidence) / 2)

    # Step 2: Calculate the Design Effect (DEFF) from network properties
    # Correction factor due to scale-free degree distribution (γ)
    # This is valid for 2 < gamma < 3
    c_gamma = (gamma - 1) / (gamma - 2)

    # Correction factor due to clustered neighborhood similarity (α)
    # This is valid for 2 < alpha < 3
    c_alpha = (alpha - 1) / (alpha - 2)

    # The total Design Effect is the product of the individual factors
    deff = c_gamma * c_alpha

    # Step 3: Formulate and calculate the sampling ratio r
    # The formula relates r to the DEFF and the statistical requirements
    # r = DEFF * (epsilon / Z)^2
    sampling_ratio = deff * (epsilon / z_score)**2

    # Step 4: Round the result to 4 decimal places
    r_rounded = round(sampling_ratio, 4)

    # Output the explanation and the final result
    print("Step-by-step calculation:")
    print(f"1. Confidence Level = {confidence} -> Z-score = {z_score:.4f}")
    print(f"2. Power-law exponent γ = {gamma} -> Correction Factor C_γ = ({gamma} - 1) / ({gamma} - 2) = {c_gamma:.4f}")
    print(f"3. Pareto shape α = {alpha} -> Correction Factor C_α = ({alpha} - 1) / ({alpha} - 2) = {c_alpha:.4f}")
    print(f"4. Total Design Effect (DEFF) = C_γ * C_α = {c_gamma:.4f} * {c_alpha:.4f} = {deff:.4f}")
    print(f"5. Tolerance ε = {epsilon}")
    print("\nFinal Equation for the sampling ratio r:")
    print(f"r = DEFF * (ε / Z)²")
    print(f"r = {deff:.4f} * ({epsilon} / {z_score:.4f})²")
    print(f"r ≈ {sampling_ratio:.8f}")
    print(f"\nThe minimum required sampling ratio r, rounded to 4 decimal places, is:")
    print(r_rounded)

    return r_rounded

# Run the calculation
final_r = calculate_sampling_ratio()
print(f"\n<<<r = {final_r}>>>")
