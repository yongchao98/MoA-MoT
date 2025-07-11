import math

def calculate_particle_size():
    """
    Calculates the particle size using the Scherrer equation based on XRD data.
    """
    # --- Given and Estimated Parameters ---
    # Scherrer constant (dimensionless)
    K = 0.9
    # X-ray wavelength in nanometers (nm)
    lambda_val = 0.1541
    # Peak position from XRD graph in degrees (2-theta) for the (111) peak
    two_theta_deg = 38.7
    # Estimated Full Width at Half Maximum (FWHM) from the graph in degrees
    beta_deg = 0.6

    # --- Calculations ---
    # Convert angles from degrees to radians for calculation
    # Bragg angle (theta) is half of the 2-theta value
    theta_rad = math.radians(two_theta_deg / 2)
    # FWHM (beta) in radians
    beta_rad = math.radians(beta_deg)

    # Bragg's angle (theta) in degrees
    theta_deg = two_theta_deg / 2

    # Scherrer equation: D = (K * lambda) / (beta * cos(theta))
    particle_size_nm = (K * lambda_val) / (beta_rad * math.cos(theta_rad))

    # --- Output Results ---
    print("--- Particle Size Calculation using Scherrer Equation ---")
    print(f"Scherrer Constant (K): {K}")
    print(f"X-ray Wavelength (λ): {lambda_val} nm")
    print(f"Peak Position (2θ): {two_theta_deg} degrees")
    print(f"Bragg Angle (θ): {theta_deg:.2f} degrees")
    print(f"Full Width at Half Maximum, FWHM (β): {beta_deg} degrees")
    print("\n--- Final Calculation ---")
    print(f"Particle Size (D) = ({K} * {lambda_val}) / ({beta_rad:.4f} rad * cos({theta_rad:.4f} rad))")
    print(f"Estimated Particle Size (D): {particle_size_nm:.2f} nm")

calculate_particle_size()
