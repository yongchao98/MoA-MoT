import math

def calculate_it3_gaussian_vertical_slice(x: float, u: float, c: float, sigma_x: float, sigma_v: float) -> float:
    """
    Calculates the upper bound of a vertical cross-section for a Gaussian Interval Type-3 Membership Function.

    This function represents the mathematical formulation that characterizes the vertical
    cross-sections by modeling the uncertainty bounds using Gaussian paradigms.

    Args:
        x: The value of the primary input variable.
        u: The value of the secondary variable (primary membership grade).
        c: The center of the Gaussian function along the x-axis.
        sigma_x: The standard deviation of the Gaussian for the reference FOU upper bound.
        sigma_v: The standard deviation of the Gaussian for the vertical slice uncertainty.

    Returns:
        The value of the upper bound of the vertical cross-section, μ̄_Ã(x, u).
    """
    # --- Step 1: Calculate the upper bound of the reference Type-2 FOU ---
    # This is the mean of the Gaussian function for the vertical slice.
    # Formula: μ̄_FOU(x) = exp(-0.5 * ((x - c) / σ_x)²)
    try:
        mu_fou_upper = math.exp(-0.5 * ((x - c) / sigma_x)**2)
    except ZeroDivisionError:
        print("Error: Standard deviation sigma_x cannot be zero.")
        return 0.0

    # --- Step 2: Calculate the upper bound of the IT3 MF vertical cross-section ---
    # This models the uncertainty in the vertical dimension around the reference FOU bound.
    # Formula: μ̄_Ã(x, u) = exp(-0.5 * ((u - μ̄_FOU(x)) / σ_v)²)
    try:
        mu_A_upper_xu = math.exp(-0.5 * ((u - mu_fou_upper) / sigma_v)**2)
    except ZeroDivisionError:
        print("Error: Standard deviation sigma_v cannot be zero.")
        return 0.0

    # --- Output the explanation and results ---
    print("--- Mathematical Formulation of the Vertical Cross-Section ---")
    print("The upper bound of the vertical cross-section, μ̄_Ã(x, u), is defined as:")
    print("μ̄_Ã(x, u) = exp( -0.5 * ( (u - μ̄_FOU(x)) / σ_v )² )\n")
    print("Where the reference FOU upper bound, μ̄_FOU(x), is:")
    print("μ̄_FOU(x) = exp( -0.5 * ( (x - c) / σ_x )² )\n")

    print("--- Calculation with Given Parameters ---")
    print(f"Parameters:")
    print(f"  Primary Input (x)         = {x}")
    print(f"  Secondary Input (u)       = {u}")
    print(f"  Center (c)                = {c}")
    print(f"  Std. Dev. for FOU (σ_x)   = {sigma_x}")
    print(f"  Std. Dev. Vertical (σ_v)  = {sigma_v}\n")

    print("Step 1: Calculate the reference FOU upper bound μ̄_FOU(x)")
    print(f"μ̄_FOU({x}) = exp( -0.5 * ( ({x} - {c}) / {sigma_x} )² )")
    print(f"μ̄_FOU({x}) = {mu_fou_upper:.4f}\n")

    print("Step 2: Calculate the IT3 vertical cross-section upper bound μ̄_Ã(x, u)")
    print(f"μ̄_Ã({x}, {u}) = exp( -0.5 * ( ({u} - {mu_fou_upper:.4f}) / {sigma_v} )² )")
    print(f"μ̄_Ã({x}, {u}) = {mu_A_upper_xu:.4f}\n")
    
    print("--- Final Result ---")
    print(f"The resulting value for the upper bound of the vertical cross-section is: {mu_A_upper_xu:.4f}")

    return mu_A_upper_xu

if __name__ == '__main__':
    # --- Define Parameters for the IT3 Membership Function ---
    # You can modify these values to explore the formulation.
    primary_input_x = 2.5
    secondary_input_u = 0.8
    gaussian_center_c = 5.0
    stdev_x = 2.0  # Std. deviation for the primary membership function shape
    stdev_v = 0.1  # Std. deviation for the vertical slice (uncertainty in u)

    # --- Execute the calculation ---
    final_value = calculate_it3_gaussian_vertical_slice(
        x=primary_input_x,
        u=secondary_input_u,
        c=gaussian_center_c,
        sigma_x=stdev_x,
        sigma_v=stdev_v
    )
