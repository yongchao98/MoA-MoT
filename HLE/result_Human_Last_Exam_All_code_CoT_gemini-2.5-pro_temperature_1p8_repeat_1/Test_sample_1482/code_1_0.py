import math

def calculate_it3_mf_vertical_cross_section():
    """
    Calculates and displays the formulation for a vertical cross-section of a 
    Gaussian Interval Type-3 Membership Function's upper bound.

    This function demonstrates the mathematical relationship between the primary
    input (x), the secondary input (u), and the resulting membership grade.
    """
    # --- Parameters ---
    # x: The fixed value of the primary input variable for the vertical cross-section.
    x_prime = 5.0
    # u: A specific value for the secondary variable, u is in [0, 1].
    # It scales the amplitude of the Gaussian function.
    u = 0.8
    # c: The center (mean) of the Gaussian function.
    c = 6.0
    # sigma_upper: The standard deviation for the upper bound of the membership function.
    # A larger sigma corresponds to greater uncertainty.
    sigma_upper = 2.0
    
    # --- Calculation ---
    # The vertical cross-section formulation for the upper bound is:
    # µ_upper(x, u) = u * exp(-0.5 * ((x - c) / sigma_upper)^2)
    
    exponent = -0.5 * math.pow((x_prime - c) / sigma_upper, 2)
    gaussian_val = math.exp(exponent)
    membership_value = u * gaussian_val
    
    # --- Output ---
    # Print the full equation with the values substituted in.
    print("The mathematical formulation for the upper bound of the vertical cross-section is:")
    print(f"µ_upper(x, u) = u * exp(-0.5 * ((x - c) / σ_upper)²) \n")
    print("Substituting the given values:")
    print(f"µ_upper({x_prime}, {u}) = {u} * exp(-0.5 * (({x_prime} - {c}) / {sigma_upper})²)")
    print(f"µ_upper({x_prime}, {u}) = {u} * exp({exponent:.4f})")
    print(f"µ_upper({x_prime}, {u}) = {u} * {gaussian_val:.4f}")
    print(f"µ_upper({x_prime}, {u}) = {membership_value:.4f}")

# Execute the function to show the result.
calculate_it3_mf_vertical_cross_section()
<<<µ_upper(5.0, 0.8) = 0.7047>>>