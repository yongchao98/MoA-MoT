import math

def calculate_vertical_slice_formulation(x_prime):
    """
    Calculates and prints the mathematical formulation for the vertical
    cross-sections of a Gaussian-based IT3 MF.

    The vertical cross-section at a fixed x=x_prime is a Type-2 Fuzzy Set
    over a secondary variable u. This T2 FS is defined by its Footprint of
    Uncertainty (FOU), bounded by a lower and an upper membership function.
    These functions are modeled as Gaussians over u, where their parameters
    (height, center, sigma) are functions of x_prime.
    """

    # --- Parameters for the UPPER bound formulation ---
    # Height of the Gaussian is itself a Gaussian function of x
    upper_height_center = 4.0
    upper_height_sigma = 2.0
    # Center of the Gaussian is a linear function of x: a*x + b
    upper_center_a = 0.05
    upper_center_b = 0.4
    # Standard deviation of the Gaussian is a linear function of x: a*x + b
    upper_sigma_a = -0.01
    upper_sigma_b = 0.2

    # --- Parameters for the LOWER bound formulation ---
    # Height of the Gaussian is a scaled Gaussian function of x
    lower_height_scale = 0.8
    lower_height_center = 6.0
    lower_height_sigma = 3.0
    # Center of the Gaussian is a linear function of x: a*x + b
    lower_center_a = 0.04
    lower_center_b = 0.3
    # Standard deviation of the Gaussian is a linear function of x: a*x + b
    lower_sigma_a = -0.005
    lower_sigma_b = 0.15

    # --- Calculation for the UPPER bound at x_prime ---
    h_upper = math.exp(-0.5 * math.pow((x_prime - upper_height_center) / upper_height_sigma, 2))
    # Ensure center for u is within [0, 1]
    c_upper = max(0, min(1, upper_center_a * x_prime + upper_center_b))
    # Ensure sigma is positive
    s_upper = max(1e-6, upper_sigma_a * x_prime + upper_sigma_b)

    # --- Calculation for the LOWER bound at x_prime ---
    h_lower_unscaled = math.exp(-0.5 * math.pow((x_prime - lower_height_center) / lower_height_sigma, 2))
    h_lower = lower_height_scale * h_lower_unscaled
    # Ensure center for u is within [0, 1]
    c_lower = max(0, min(1, lower_center_a * x_prime + lower_center_b))
    # Ensure sigma is positive
    s_lower = max(1e-6, lower_sigma_a * x_prime + lower_sigma_b)


    # --- Print the final formulations ---
    print("For a fixed primary input x' = {:.1f}, the vertical cross-section is a Type-2 fuzzy set over the secondary variable u.".format(x_prime))
    print("Its Footprint of Uncertainty (FOU) is bounded by the following Gaussian-based membership functions:\n")

    # Final Upper Bound Equation
    # Using python's string formatting to inject the calculated numbers.
    # Each number is explicitly included as requested.
    upper_equation = "μ_upper(u) = {:.3f} * exp(-0.5 * ((u - {:.3f}) / {:.3f})**2)".format(h_upper, c_upper, s_upper)
    print("Upper Bound Formulation:")
    print(upper_equation)
    print("\n")

    # Final Lower Bound Equation
    # Each number is explicitly included as requested.
    lower_equation = "μ_lower(u) = {:.3f} * exp(-0.5 * ((u - {:.3f}) / {:.3f})**2)".format(h_lower, c_lower, s_lower)
    print("Lower Bound Formulation:")
    print(lower_equation)


if __name__ == '__main__':
    # Define the fixed value for the primary input variable x
    fixed_x = 5.0
    calculate_vertical_slice_formulation(fixed_x)
