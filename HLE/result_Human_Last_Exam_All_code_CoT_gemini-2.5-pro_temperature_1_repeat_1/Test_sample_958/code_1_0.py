import numpy as np

def calculate_stored_energy(n, theta_deg, omega_over_c=1.0, E_x0_i_sq=1.0):
    """
    Calculates the time-averaged stored energy per unit area for an evanescent wave
    based on the formulas in Option D.

    Args:
        n (float): Refractive index of the core material.
        theta_deg (float): Incident angle in degrees.
        omega_over_c (float, optional): Ratio of angular frequency to the speed of light. Defaults to 1.0.
        E_x0_i_sq (float, optional): The square of the magnitude of the x-component of the
                                     incident electric field amplitude. Defaults to 1.0.
    """
    # Physical constants
    eps0 = 8.854e-12  # Permittivity of free space in F/m

    # Convert angle to radians
    theta_rad = np.deg2rad(theta_deg)
    
    # Check for total internal reflection condition
    if n * np.sin(theta_rad) <= 1:
        print("Total internal reflection condition not met. Evanescent wave is not generated.")
        return

    # Intermediate terms for clarity
    s2 = np.sin(theta_rad)**2
    n2 = n**2
    
    # Numerators from Option D
    numerator_E = n2 * (2 * n2 * s2 - 1)
    numerator_H = n2 * (n2 * s2 - 1)

    # Denominator (common for both E and H field energy)
    sqrt_term = np.sqrt(n2 * s2 - 1)
    denom_factor1 = 2 * omega_over_c
    denom_factor2 = (n2 - 1)
    denom_factor3 = ((n2 + 1) * s2 - 1)
    denominator = denom_factor1 * denom_factor2 * denom_factor3 * sqrt_term
    
    # Check if denominator is zero
    if denominator == 0:
        print("Denominator is zero, cannot calculate energy.")
        return

    # Calculate stored energy per unit area
    energy_E_factor = numerator_E / denominator
    energy_H_factor = numerator_H / denominator
    
    energy_E = energy_E_factor * eps0 * E_x0_i_sq
    energy_H = energy_H_factor * eps0 * E_x0_i_sq

    print("--- Calculation based on Option D ---")
    print(f"Inputs: n = {n}, theta = {theta_deg} degrees")
    print("\n--- Equation for Energy in E field ---")
    print(f"Energy_E = (n^2 * (2*n^2*sin^2(θ) - 1)) / (2*(ω/c)*(n^2-1)*((n^2+1)*sin^2(θ)-1)*sqrt(n^2*sin^2(θ)-1)) * ε₀*|E_x0_i|²")
    print(f"Numerator term: n^2 * (2*n^2*sin^2(θ) - 1) = {numerator_E:.4f}")
    print(f"Denominator term: 2*(ω/c)*(n^2-1)*((n^2+1)*sin^2(θ)-1)*sqrt(n^2*sin^2(θ)-1) = {denominator:.4f}")
    print(f"Resulting Energy in E-field = {energy_E_factor:.4f} * ε₀ * |E_x0_i|² = {energy_E:.4e} J/m² (for |E_x0_i|²=1, ω/c=1)")
    
    print("\n--- Equation for Energy in H field ---")
    print(f"Energy_H = (n^2 * (n^2*sin^2(θ) - 1)) / (2*(ω/c)*(n^2-1)*((n^2+1)*sin^2(θ)-1)*sqrt(n^2*sin^2(θ)-1)) * ε₀*|E_x0_i|²")
    print(f"Numerator term: n^2 * (n^2*sin^2(θ) - 1) = {numerator_H:.4f}")
    print(f"Denominator term: 2*(ω/c)*(n^2-1)*((n^2+1)*sin^2(θ)-1)*sqrt(n^2*sin^2(θ)-1) = {denominator:.4f}")
    print(f"Resulting Energy in H-field = {energy_H_factor:.4f} * ε₀ * |E_x0_i|² = {energy_H:.4e} J/m² (for |E_x0_i|²=1, ω/c=1)")


# Example usage:
# Glass core (n=1.5) and an incident angle of 60 degrees.
# The critical angle is asin(1/1.5) = 41.8 degrees, so 60 degrees causes TIR.
calculate_stored_energy(n=1.5, theta_deg=60)