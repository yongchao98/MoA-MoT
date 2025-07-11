import numpy as np

def calculate_evanescent_wave_energy():
    """
    Calculates the time-averaged stored energy per unit area for the electric and
    magnetic fields of an evanescent wave for a p-polarized light.

    The function uses the formulas provided in option D, which are:
    Energy in E field = (n^2*(2*n^2*sin^2(theta) - 1)) / Denominator * const
    Energy in H field = (n^2*(n^2*sin^2(theta) - 1)) / Denominator * const
    Denominator = 2*(omega/c)*(n^2-1)*((n^2+1)*sin^2(theta)-1)*sqrt(n^2*sin^2(theta)-1)
    Constant = epsilon_0 * |E_x0_i|^2
    """

    # --- Input Parameters (using realistic values for an optical fiber) ---
    # Refractive index of the core
    n = 1.45
    # Angle of incidence in degrees (must be > critical angle)
    theta_deg = 60.0
    # Wavelength of light in meters (e.g., 1550 nm for telecommunications)
    wavelength = 1550e-9
    # Amplitude of the incident electric field's x-component (V/m)
    E_x0_i = 1.0  # Normalized to 1

    # --- Physical Constants ---
    # Speed of light in vacuum (m/s)
    c = 299792458.0
    # Vacuum permittivity (F/m)
    epsilon_0 = 8.854187817e-12

    # --- Derived quantities ---
    # Angle of incidence in radians
    theta = np.deg2rad(theta_deg)
    # Angular frequency (rad/s)
    omega = 2 * np.pi * c / wavelength
    # Square of sin(theta)
    s2 = np.sin(theta)**2

    # --- Check for Total Internal Reflection ---
    critical_angle_rad = np.arcsin(1/n)
    if theta <= critical_angle_rad:
        print("Warning: Angle of incidence is not greater than the critical angle.")
        print(f"Theta = {theta_deg:.2f} degrees, Critical Angle = {np.rad2deg(critical_angle_rad):.2f} degrees")
        print("Evanescent wave is not generated under these conditions.")
        return

    # --- Calculation using formulas from Option D ---
    # This factor appears in the denominator of both expressions
    sqrt_term = np.sqrt(n**2 * s2 - 1)

    # Denominator calculation
    # Denominator = 2*(omega/c)*(n^2-1)*((n^2+1)*sin^2(theta)-1)*sqrt(n^2*sin^2(theta)-1)
    term1 = 2 * (omega / c)
    term2 = n**2 - 1
    term3 = (n**2 + 1) * s2 - 1
    denominator = term1 * term2 * term3 * sqrt_term

    # Constant multiplier
    # Constant = epsilon_0 * |E_x0_i|^2
    const_multiplier = epsilon_0 * E_x0_i**2

    # Numerator for Electric Field Energy
    # Numerator_E = n^2*(2*n^2*sin^2(theta) - 1)
    numerator_E = n**2 * (2 * n**2 * s2 - 1)

    # Numerator for Magnetic Field Energy
    # Numerator_H = n^2*(n^2*sin^2(theta) - 1)
    numerator_H = n**2 * (n**2 * s2 - 1)
    
    # Final energy calculations
    energy_E = (numerator_E / denominator) * const_multiplier
    energy_H = (numerator_H / denominator) * const_multiplier

    # --- Outputting the Results ---
    print("--- Input Parameters ---")
    print(f"Refractive index of core (n): {n}")
    print(f"Angle of incidence (theta): {theta_deg} degrees")
    print(f"Wavelength (lambda): {wavelength * 1e9} nm")
    print(f"Incident E-field amplitude |E_x0_i|: {E_x0_i} V/m\n")
    
    print("--- Calculation of Stored Energy per Unit Area ---")
    print("\n--- Energy in Electric Field (W'_E) ---")
    print("Formula: (n^2 * (2*n^2*sin^2(theta) - 1)) / (2*(omega/c)*(n^2-1)*((n^2+1)*sin^2(theta)-1)*sqrt(n^2*sin^2(theta)-1)) * epsilon_0 * |E_x0_i|^2")
    print("\nNumerical values for the formula parts:")
    print(f"  n^2 = {n**2:.4f}")
    print(f"  sin^2(theta) = {s2:.4f}")
    print(f"  Numerator_E = {numerator_E:.4f}")
    print(f"  Denominator = {denominator:.4e}")
    print(f"  Constant Part (epsilon_0 * |E_x0_i|^2) = {const_multiplier:.4e}")
    print(f"\nResulting W'_E = {energy_E:.4e} J/m^2\n")

    print("\n--- Energy in Magnetic Field (W'_H) ---")
    print("Formula: (n^2 * (n^2*sin^2(theta) - 1)) / (2*(omega/c)*(n^2-1)*((n^2+1)*sin^2(theta)-1)*sqrt(n^2*sin^2(theta)-1)) * epsilon_0 * |E_x0_i|^2")
    print("\nNumerical values for the formula parts:")
    print(f"  n^2 = {n**2:.4f}")
    print(f"  sin^2(theta) = {s2:.4f}")
    print(f"  Numerator_H = {numerator_H:.4f}")
    print(f"  Denominator = {denominator:.4e}")
    print(f"  Constant Part (epsilon_0 * |E_x0_i|^2) = {const_multiplier:.4e}")
    print(f"\nResulting W'_H = {energy_H:.4e} J/m^2")

if __name__ == '__main__':
    calculate_evanescent_wave_energy()
