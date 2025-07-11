import numpy as np

def calculate_evanescent_energy(n, theta_deg, lambda_nm, E_x0_i):
    """
    Calculates the time-averaged stored energy per unit area in the evanescent wave.

    Args:
        n (float): Refractive index of the core material.
        theta_deg (float): Incident angle in degrees.
        lambda_nm (float): Wavelength of light in nanometers.
        E_x0_i (float): Amplitude of the x-component of the incident electric field.
    """
    # --- Constants ---
    # Vacuum permittivity in F/m
    epsilon_0 = 8.854187817e-12
    # Speed of light in m/s
    c = 299792458.0

    # --- Input conversions ---
    theta_rad = np.deg2rad(theta_deg)
    # Convert wavelength from nm to m
    lambda_m = lambda_nm * 1e-9
    # omega/c = 2*pi/lambda
    w_over_c = 2 * np.pi / lambda_m
    E_x0_i_sq = E_x0_i**2

    # --- Check for TIR condition ---
    if n * np.sin(theta_rad) <= 1:
        print("Total Internal Reflection condition not met (n*sin(theta) <= 1).")
        print(f"Critical angle is {np.rad2deg(np.arcsin(1/n)):.2f} degrees.")
        return

    # --- Intermediate Calculations based on Option D ---
    sin_theta_sq = np.sin(theta_rad)**2
    # This is gamma^2 in the derivation
    gamma_sq = n**2 * sin_theta_sq - 1
    # This is gamma
    sqrt_gamma_sq = np.sqrt(gamma_sq)

    # Denominator of the energy expressions
    denominator = 2 * w_over_c * (n**2 - 1) * ((n**2 + 1) * sin_theta_sq - 1) * sqrt_gamma_sq

    # --- Energy in Electric Field (W_E) ---
    numerator_E = n**2 * (2 * n**2 * sin_theta_sq - 1) * epsilon_0 * E_x0_i_sq
    W_E = numerator_E / denominator

    # --- Energy in Magnetic Field (W_H) ---
    # Note: The expression from choice D is used.
    numerator_H = n**2 * gamma_sq * epsilon_0 * E_x0_i_sq
    W_H = numerator_H / denominator

    # --- Output Results ---
    print("Assuming the formulas from Choice D are correct:")
    print(f"For n = {n}, theta = {theta_deg} deg, lambda = {lambda_nm} nm, |E_x0_i| = {E_x0_i} V/m:")
    print(f"Energy in E field = {W_E:.4g} J/m^2")
    print(f"Energy in H field = {W_H:.4g} J/m^2")

# --- Example Usage ---
# Realistic parameters for a silica fiber core and a common telecom wavelength
n_core = 1.45
# Angle must be greater than the critical angle arcsin(1/1.45) = 43.6 deg
incident_angle = 60.0
wavelength = 1550.0
# Assume a normalized incident field amplitude
E_field_amplitude = 1.0

calculate_evanescent_energy(n_core, incident_angle, wavelength, E_field_amplitude)