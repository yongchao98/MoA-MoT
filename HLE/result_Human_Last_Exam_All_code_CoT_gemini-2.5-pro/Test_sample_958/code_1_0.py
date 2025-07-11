import numpy as np

def calculate_energy_coefficients(n, theta_deg):
    """
    Calculates the numerical coefficients for the time-averaged stored energy
    per unit area in the electric and magnetic fields of an evanescent wave,
    based on the formulas in Option D.

    The energy is given by:
    Energy = Coeff * (1/(w/c)) * const * |E_x0^i|^2
    This function calculates 'Coeff'.

    Args:
        n (float): Refractive index of the core.
        theta_deg (float): Incident angle in degrees.
    """
    # Check for TIR condition
    if n <= 1:
        print("Refractive index n must be > 1.")
        return
    theta_c_rad = np.arcsin(1/n)
    theta_c_deg = np.rad2deg(theta_c_rad)
    if theta_deg <= theta_c_deg:
        print(f"Angle {theta_deg}° is not greater than the critical angle {theta_c_deg:.2f}°. No evanescent wave is generated.")
        return

    theta_rad = np.deg2rad(theta_deg)
    sin_theta_sq = np.sin(theta_rad)**2
    n_sq = n**2

    # Common terms in the denominator
    term_sqrt = np.sqrt(n_sq * sin_theta_sq - 1)
    denom_common = 2 * (n_sq - 1) * ((n_sq + 1) * sin_theta_sq - 1) * term_sqrt

    # Numerator for Electric field energy
    num_E = n_sq * (2 * n_sq * sin_theta_sq - 1)

    # Numerator for Magnetic field energy (from Option D)
    num_H = n_sq * (n_sq * sin_theta_sq - 1)
    
    # Coefficients (the numerical part of the formula, excluding (w/c), constants, and E_field^2)
    coeff_E = num_E / denom_common
    coeff_H = num_H / denom_common

    print(f"For n = {n} and theta = {theta_deg}°:")
    print("-" * 30)
    
    # Print the full formula for Energy in E field
    print("Energy in E field =")
    print(f"  {num_E:.4f} / (2 * ({n_sq - 1:.4f}) * (({n_sq + 1:.4f}) * {sin_theta_sq:.4f} - 1) * sqrt({n_sq * sin_theta_sq - 1:.4f})) * (1/(w/c)) * epsilon_0 * |E_x0^i|^2")
    print(f"= {coeff_E:.4f} * (1/(w/c)) * epsilon_0 * |E_x0^i|^2")
    print("-" * 30)

    # Print the full formula for Energy in H field
    print("Energy in H field =")
    print(f"  {num_H:.4f} / (2 * ({n_sq - 1:.4f}) * (({n_sq + 1:.4f}) * {sin_theta_sq:.4f} - 1) * sqrt({n_sq * sin_theta_sq - 1:.4f})) * (1/(w/c)) * epsilon_0 * |E_x0^i|^2")
    print(f"= {coeff_H:.4f} * (1/(w/c)) * epsilon_0 * |E_x0^i|^2")
    print("-" * 30)


# --- User Input ---
# Example values for a glass core
n_core = 1.5
# Angle greater than the critical angle (arcsin(1/1.5) = 41.8 deg)
incident_angle = 60.0

calculate_energy_coefficients(n_core, incident_angle)