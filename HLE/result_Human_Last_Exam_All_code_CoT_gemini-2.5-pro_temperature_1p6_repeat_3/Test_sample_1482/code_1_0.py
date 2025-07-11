import math

def generate_vertical_cross_section_formulation():
    """
    Calculates and prints the mathematical formulation for the upper bound of a vertical 
    cross-section of a Gaussian Interval Type-3 Membership Function (IT3 MF).
    """
    # --- Step 1: Define parameters for the IT3 MF ---
    
    # Parameters for the primary Gaussian membership functions
    primary_mean_c = 5.0
    primary_sigma_upper = 2.0  # Std. dev. for the upper membership function (UMF)
    primary_sigma_lower = 1.0  # Std. dev. for the lower membership function (LMF)
    
    # The specific primary input value for the vertical cross-section
    x_val = 4.0
    
    # Proportionality constant for the secondary uncertainty
    k = 0.5

    # --- Step 2: Calculate primary membership grades at x_val ---
    
    # Define a helper function for Gaussian calculation
    def gaussian(x, mean, sigma):
        if sigma == 0:
            return 1.0 if x == mean else 0.0
        exponent = -((x - mean)**2) / (2 * sigma**2)
        return math.exp(exponent)

    mu_upper_at_x = gaussian(x_val, primary_mean_c, primary_sigma_upper)
    mu_lower_at_x = gaussian(x_val, primary_mean_c, primary_sigma_lower)
    
    # --- Step 3: Calculate parameters for the vertical cross-section's UMF ---
    
    # The mean of the secondary Gaussian is set to the primary UMF value
    secondary_mean_u = mu_upper_at_x
    
    # The std. dev. of the secondary Gaussian is proportional to the primary FOU height
    secondary_std_dev_u = k * (mu_upper_at_x - mu_lower_at_x)

    # --- Step 4: Construct and print the final mathematical formulation ---
    
    # We formulate the upper bound of the vertical cross-section at x_val=4.0
    # The formula is: f_x(u) = exp( -(u - secondary_mean_u)^2 / (2 * secondary_std_dev_u^2) )
    
    print("--- Interval Type-3 Fuzzy Logic System: Vertical Cross-Section Formulation ---")
    print(f"\nGiven primary input x = {x_val}")
    print(f"Primary UMF at x: overline_mu({x_val}) = {mu_upper_at_x:.4f}")
    print(f"Primary LMF at x: underline_mu({x_val}) = {mu_lower_at_x:.4f}")
    print(f"\nThe derived parameters for the vertical cross-section's UMF are:")
    print(f"  - Mean (center) of the secondary Gaussian: c_u = {secondary_mean_u:.4f}")
    print(f"  - Std. Deviation of the secondary Gaussian: sigma_u = {secondary_std_dev_u:.4f}")
    
    secondary_std_dev_u_sq = secondary_std_dev_u**2
    denominator = 2 * secondary_std_dev_u_sq

    # Final Equation string
    final_equation = f"overline_f_x={x_val}(u) = exp( -(u - {secondary_mean_u:.4f})^2 / (2 * {secondary_std_dev_u:.4f}^2) )"
    final_equation_simplified = f"overline_f_x={x_val}(u) = exp( -(u - {secondary_mean_u:.4f})^2 / {denominator:.4f} )"

    print("\nThe mathematical formulation for the Upper Membership Function of the vertical cross-section is:")
    print(final_equation)
    print("\nSimplified with denominator calculated:")
    print(final_equation_simplified)

    return final_equation_simplified

# Execute the function to get the result
final_formula = generate_vertical_cross_section_formulation()

# The final answer format as requested
# The content is the final simplified equation.
# final_formula = <<<...>>>