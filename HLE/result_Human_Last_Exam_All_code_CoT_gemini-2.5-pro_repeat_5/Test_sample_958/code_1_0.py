import numpy as np

def calculate_evanescent_wave_energy(n, theta_deg, omega, E_x0_i):
    """
    Calculates the time-averaged stored energy per unit area for an evanescent wave.

    Args:
        n (float): Refractive index of the core material.
        theta_deg (float): Incident angle in degrees.
        omega (float): Angular frequency of the light in rad/s.
        E_x0_i (float): Amplitude of the x-component of the incident electric field in V/m.
    """
    # Physical constants
    epsilon_0 = 8.854e-12  # Permittivity of free space in F/m
    c = 2.998e8          # Speed of light in m/s

    # Convert angle to radians
    theta = np.deg2rad(theta_deg)

    # Check for TIR condition
    theta_c = np.arcsin(1 / n)
    if theta <= theta_c:
        print("Warning: Incident angle is not greater than the critical angle.")
        print(f"Critical angle is {np.rad2deg(theta_c):.2f} degrees.")
        return

    # Intermediate terms for clarity
    k0 = omega / c
    n2 = n**2
    sin_theta = np.sin(theta)
    sin2_theta = sin_theta**2
    
    term_sqrt = np.sqrt(n2 * sin2_theta - 1)
    
    # Denominator from the formula
    # Denominator = 2 * (omega/c) * (n^2 - 1) * [(n^2 + 1)sin^2(theta) - 1] * sqrt(n^2*sin^2(theta) - 1)
    denominator = 2 * k0 * (n2 - 1) * ((n2 + 1) * sin2_theta - 1) * term_sqrt

    # Numerator for Electric Field Energy
    # Num_E = n^2 * (2*n^2*sin^2(theta) - 1)
    numerator_E = n2 * (2 * n2 * sin2_theta - 1)

    # Numerator for Magnetic Field Energy
    # Num_H = n^2 * (n^2*sin^2(theta) - 1)
    numerator_H = n2 * (n2 * sin2_theta - 1)

    # Final calculation
    common_factor = epsilon_0 * (E_x0_i**2) / denominator
    
    energy_E = numerator_E * common_factor
    energy_H = numerator_H * common_factor

    print("--- Calculated Stored Energy per Unit Area ---")
    print(f"For n = {n}, theta = {theta_deg} deg:")
    print(f"Energy in E field = {energy_E:.4e} J/m^2")
    print(f"Energy in H field = {energy_H:.4e} J/m^2")
    print("\n--- Based on Formula from choice D ---")
    print("Energy in E field = (n^2 * (2*n^2*sin^2(theta) - 1)) / Denominator * epsilon_0 * |E_x0^i|^2")
    print("Energy in H field = (n^2 * (n^2*sin^2(theta) - 1)) / Denominator * epsilon_0 * |E_x0^i|^2")
    print("where Denominator = 2*(omega/c)*(n^2-1)*((n^2+1)*sin^2(theta)-1)*sqrt(n^2*sin^2(theta)-1)")


# Example Usage:
# Parameters for a typical silica fiber and HeNe laser
n_core = 1.458  # Refractive index of silica core
wavelength = 632.8e-9  # Wavelength of HeNe laser in meters
angular_freq = 2 * np.pi * c / wavelength
incident_angle = 68.0  # Angle in degrees (critical angle is ~43.3 deg)
E_field_amplitude = 1.0  # Normalized incident field amplitude in V/m

calculate_evanescent_wave_energy(n_core, incident_angle, angular_freq, E_field_amplitude)
