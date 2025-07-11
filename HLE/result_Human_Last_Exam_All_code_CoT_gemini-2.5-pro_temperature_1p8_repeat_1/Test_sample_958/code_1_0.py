import numpy as np

def calculate_evanescent_energy(n, theta_deg):
    """
    Calculates the prefactors for the time-averaged stored energy per unit area
    in an evanescent wave for p-polarized light.

    The calculation is based on the formulas:
    Energy_E = C_E * (epsilon_0 * |E_x0_i|^2) / (omega/c)
    Energy_H = C_H * (epsilon_0 * |E_x0_i|^2) / (omega/c)
    (Note: The problem's formula for H_energy has epsilon_0, not mu_0)

    Args:
        n (float): Refractive index of the core material.
        theta_deg (float): Incident angle in degrees.

    Returns:
        tuple: A tuple containing the dimensionless prefactors for electric (C_E)
               and magnetic (C_H) energy. Returns (None, None) if TIR condition
               is not met.
    """
    # Critical angle for total internal reflection (TIR)
    if n <= 1:
        print("Error: Refractive index n must be > 1.")
        return None, None
    theta_c_deg = np.degrees(np.arcsin(1 / n))

    if theta_deg <= theta_c_deg:
        print(f"Angle {theta_deg} deg is not greater than the critical angle {theta_c_deg:.2f} deg.")
        print("No evanescent wave is generated under these conditions.")
        return None, None

    # Convert angle to radians for trigonometric functions
    theta_rad = np.radians(theta_deg)

    # Pre-calculate common terms
    n2 = n**2
    sin_theta = np.sin(theta_rad)
    sin2_theta = sin_theta**2

    # Check for valid denominator
    common_denom_sqrt_term = n2 * sin2_theta - 1
    if common_denom_sqrt_term <= 0:
        # This case is already handled by the critical angle check but is good for robustness
        return None, None
    sqrt_term = np.sqrt(common_denom_sqrt_term)

    common_denominator = 2 * (n2 - 1) * ((n2 + 1) * sin2_theta - 1) * sqrt_term

    if common_denominator == 0:
        print("Error: Denominator is zero.")
        return None, None

    # Numerator for the Electric Field Energy prefactor
    num_E = n2 * (2 * n2 * sin2_theta - 1)

    # Numerator for the Magnetic Field Energy prefactor
    num_H = n2 * (n2 * sin2_theta - 1)

    # Prefactors
    C_E = num_E / common_denominator
    C_H = num_H / common_denominator
    
    # Print the formulas and the results
    print("Based on the derived formulas corresponding to choice D:")
    print("\nEnergy in E field = [n^2 * (2*n^2*sin^2(θ) - 1)] / [2*(n^2-1)*((n^2+1)*sin^2(θ)-1)*sqrt(n^2*sin^2(θ)-1)] * (ε₀*|E_x₀ⁱ|²) / (ω/c)")
    print(f"Energy in E field = [{n2:.2f} * (2*{n2:.2f}*{sin2_theta:.3f} - 1)] / [2*({n2-1:.2f})*(({n2+1:.2f})*{sin2_theta:.3f}-1)*{sqrt_term:.3f}] * (ε₀*|E_x₀ⁱ|²) / (ω/c)")
    print(f"Calculated Prefactor C_E: {C_E:.4f}\n")
    
    print("Energy in H field = [n^2 * (n^2*sin^2(θ) - 1)] / [2*(n^2-1)*((n^2+1)*sin^2(θ)-1)*sqrt(n^2*sin^2(θ)-1)] * (ε₀*|E_x₀ⁱ|²))")
    print(f"Energy in H field = [{n2:.2f} * ({n2:.2f}*{sin2_theta:.3f} - 1)] / [2*({n2-1:.2f})*(({n2+1:.2f})*{sin2_theta:.3f}-1)*{sqrt_term:.3f}] * (ε₀*|E_x₀ⁱ|²)")
    print(f"Calculated Prefactor C_H: {C_H:.4f}\n")


# --- Main execution ---
# Example parameters for a typical glass core fiber
n_core = 1.5
theta_incident_deg = 60.0

calculate_evanescent_energy(n_core, theta_incident_deg)
