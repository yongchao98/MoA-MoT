import math

def generate_vertical_slice_formulation():
    """
    Calculates and prints the mathematical formulation for the vertical cross-section
    of a Gaussian-based Interval Type-3 Membership Function (IT3 MF).
    """

    # --- Step 1: Define parameters for the Gaussian functions ---
    # Parameters for the primary Gaussian function in the 'x' domain
    c_x = 5.0      # Mean of the primary Gaussian
    sigma_x = 2.0  # Standard deviation of the primary Gaussian

    # Parameter for the secondary Gaussian function in the 'u' domain
    sigma_u = 0.1  # Standard deviation of the secondary Gaussian, controlling uncertainty

    # The fixed primary input variable for the vertical cross-section
    x_prime = 4.0

    # --- Step 2: Calculate the mean for the secondary Gaussian ---
    # The mean of the secondary Gaussian is the value of the primary Gaussian at x = x'
    exponent_primary = -((x_prime - c_x)**2) / (2 * sigma_x**2)
    mean_secondary = math.exp(exponent_primary)

    # --- Step 3: Print the formulation with calculated numbers ---
    print("This script determines the mathematical formulation for the upper bound of a vertical cross-section of a Gaussian IT3 MF.")
    print("-" * 80)
    print("The general formulation for the upper bound of the cross-section at x' is:")
    print("  overline_μ(u | x=x') = exp( -(u - μ_P(x'))^2 / (2 * σ_u^2) )")
    print("\nWhere μ_P(x') is the value of the primary Gaussian at x'.")
    print("-" * 80)
    print("\nGiven the following parameters:")
    print(f"  - Fixed primary input, x' = {x_prime}")
    print(f"  - Primary Gaussian mean, c = {c_x}")
    print(f"  - Primary Gaussian std dev, σ_x = {sigma_x}")
    print(f"  - Secondary Gaussian std dev, σ_u = {sigma_u}")

    print("\nStep-by-step calculation:")
    print(f"  1. The mean of the secondary Gaussian is μ_P({x_prime}) = exp( -({x_prime} - {c_x})^2 / (2 * {sigma_x}^2) ) = {mean_secondary:.4f}")
    
    print("\nFinal Equation for the Vertical Cross-Section's Upper Bound:")
    print("The final formulation, with all numbers substituted, is:")
    # This print statement shows each individual number in the equation.
    print(f"  overline_μ(u) = exp( -(u - {mean_secondary:.4f})^2 / (2 * {sigma_u}^2) )")

generate_vertical_slice_formulation()

<<<overline_μ(u) = exp( -(u - 0.8825)^2 / (2 * 0.1^2) )>>>