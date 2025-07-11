import numpy as np

def calculate_evanescent_wave_energy(n, theta_deg, omega_c, E_xi_sq, epsilon_0=8.854e-12):
    """
    Calculates the time-averaged stored energy per unit area for an evanescent wave.

    Args:
        n (float): Refractive index of the core material.
        theta_deg (float): Incident angle in degrees.
        omega_c (float): Angular frequency divided by the speed of light (omega/c).
        E_xi_sq (float): The squared magnitude of the x-component of the incident electric field, |E_x0^i|^2.
        epsilon_0 (float): Permittivity of free space.

    Returns:
        tuple: A tuple containing the energy in the E field and H field.
    """
    # Convert angle to radians
    theta_rad = np.deg2rad(theta_deg)
    s = np.sin(theta_rad)
    s2 = s**2
    n2 = n**2

    # Check for TIR condition
    sin_theta_c = 1 / n
    if s <= sin_theta_c:
        print("Warning: Angle is not greater than the critical angle. No evanescent wave is formed.")
        return 0, 0
    
    # Calculate terms in the formula from choice D
    sqrt_term = np.sqrt(n2 * s2 - 1)
    
    # Common denominator part of the expressions
    # Note: the 2*omega_c is part of the expression
    common_denominator = (2 * omega_c * (n2 - 1) * 
                          ((n2 + 1) * s2 - 1) * sqrt_term)

    if common_denominator == 0:
        return float('inf'), float('inf')

    # Numerator for Electric Field Energy
    numerator_E = n2 * (2 * n2 * s2 - 1)
    
    # Numerator for Magnetic Field Energy
    numerator_H = n2 * (n2 * s2 - 1)
    
    # Calculate energies
    energy_E = (numerator_E / common_denominator) * epsilon_0 * E_xi_sq
    energy_H = (numerator_H / common_denominator) * epsilon_0 * E_xi_sq
    
    # To display the formula structure in the output
    print("Based on the formula from the correct answer choice:")
    print(f"Energy in E field = (n^2 * (2*n^2*sin^2(theta) - 1)) / (2*(w/c)*(n^2-1)*((n^2+1)*sin^2(theta)-1)*sqrt(n^2*sin^2(theta)-1)) * epsilon_0 * |E_x0^i|^2")
    print("Plugging in the numbers:")
    print(f"Energy in E field = ({numerator_E:.4f}) / ({common_denominator:.4f}) * epsilon_0 * |E_x0^i|^2")
    print(f"                   = {energy_E:.4e} J/m^2 (assuming |E_x0^i|^2 = 1, w/c = 1)")
    print("-" * 20)
    print(f"Energy in H field = (n^2 * (n^2*sin^2(theta) - 1)) / (2*(w/c)*(n^2-1)*((n^2+1)*sin^2(theta)-1)*sqrt(n^2*sin^2(theta)-1)) * epsilon_0 * |E_x0^i|^2")
    print("Plugging in the numbers:")
    print(f"Energy in H field = ({numerator_H:.4f}) / ({common_denominator:.4f}) * epsilon_0 * |E_x0^i|^2")
    print(f"                   = {energy_H:.4e} J/m^2 (assuming |E_x0^i|^2 = 1, w/c = 1)")
    
    return energy_E, energy_H

# --- Example Usage ---
# Parameters for the fiber-optic system
n_core = 1.5           # Refractive index of the core
incident_angle = 60.0  # Incident angle in degrees > critical angle (41.8 deg)
omega_over_c = 1.0e7   # Example value for omega/c (corresponds to lambda ~ 188 nm)
E_x0_i_sq = 1.0        # Normalized squared amplitude of incident E-field x-component

print("Calculating stored energy for an example case:\n")
calculate_evanescent_wave_energy(n_core, incident_angle, omega_over_c, E_x0_i_sq)