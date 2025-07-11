import numpy as np

def calculate_evanescent_wave_energy(n, theta_deg, omega_c, E_x0_i_sq, epsilon_0):
    """
    Calculates the time-averaged stored energy per unit area for the evanescent wave.
    
    Args:
        n (float): Refractive index of the core material.
        theta_deg (float): Incident angle in degrees.
        omega_c (float): Angular frequency divided by the speed of light (omega/c).
        E_x0_i_sq (float): The squared magnitude of the incident electric field's x-component |E_x0^i|^2.
        epsilon_0 (float): Vacuum permittivity.
    """
    
    # Convert angle to radians
    theta_rad = np.deg2rad(theta_deg)
    
    # Check for Total Internal Reflection condition
    if n * np.sin(theta_rad) <= 1:
        print("Warning: Total Internal Reflection condition (n*sin(theta) > 1) is not met.")
        print(f"n*sin(theta) = {n * np.sin(theta_rad):.4f}")
        return

    # Calculate recurring terms
    sin_sq_theta = np.sin(theta_rad)**2
    
    # Expression based on Answer Choice D
    # Denominator components
    denom_term1 = np.sqrt(n**2 * sin_sq_theta - 1)
    denom_term2 = (n**2 - 1)
    denom_term3 = (n**2 + 1) * sin_sq_theta - 1
    
    # Full denominator
    denominator = 2 * omega_c * denom_term2 * denom_term3 * denom_term1

    # --- Energy in Electric Field (W_e) ---
    # Numerator for W_e
    num_e_term = 2 * n**2 * sin_sq_theta - 1
    numerator_e = n**2 * num_e_term * epsilon_0 * E_x0_i_sq
    
    # Calculate W_e
    W_e = numerator_e / denominator

    # --- Energy in Magnetic Field (W_h) ---
    # Numerator for W_h
    num_h_term = n**2 * sin_sq_theta - 1
    numerator_h = n**2 * num_h_term * epsilon_0 * E_x0_i_sq

    # Calculate W_h
    W_h = numerator_h / denominator

    # --- Print results ---
    print("Based on the formulas in Choice D:")
    print(f"Energy in E field = (n^2 * (2*n^2*sin^2(theta) - 1)) / (2 * (omega/c) * (n^2-1) * ((n^2+1)*sin^2(theta)-1) * sqrt(n^2*sin^2(theta)-1)) * epsilon_0 * |E_x0_i|^2")
    print(f"Energy in H field = (n^2 * (n^2*sin^2(theta) - 1)) / (2 * (omega/c) * (n^2-1) * ((n^2+1)*sin^2(theta)-1) * sqrt(n^2*sin^2(theta)-1)) * epsilon_0 * |E_x0_i|^2")
    print("-" * 30)
    print("Given values:")
    print(f"  n = {n}")
    print(f"  theta = {theta_deg} degrees")
    print(f"  omega/c = {omega_c:.2e} m^-1")
    print(f"  |E_x0_i|^2 = {E_x0_i_sq} (V/m)^2")
    print(f"  epsilon_0 = {epsilon_0:.4e} F/m")
    print("-" * 30)
    print("Calculated energy per unit area:")
    print(f"Energy in E field: {W_e:.4e} J/m^2")
    print(f"Energy in H field: {W_h:.4e} J/m^2")


if __name__ == '__main__':
    # Define physical constants and parameters for a numerical example
    # These values are chosen for demonstration purposes.
    n_core = 1.5           # Refractive index of the core (e.g., glass)
    theta_incident = 60.0  # Incident angle in degrees (must be > critical angle)
    
    # For n=1.5, critical angle is asin(1/1.5) = 41.8 degrees. 60 degrees is valid.
    
    # Example light parameters
    wavelength = 1550e-9   # Wavelength in meters (common for fiber optics)
    c = 3.0e8              # Speed of light in vacuum (m/s)
    omega = 2 * np.pi * c / wavelength
    omega_over_c = omega / c # This is 2*pi/wavelength
    
    E_x0_i_sq_val = 1.0    # Normalized squared amplitude of the incident E-field x-component
    epsilon_0_val = 8.854e-12 # Vacuum permittivity in F/m
    
    calculate_evanescent_wave_energy(n_core, theta_incident, omega_over_c, E_x0_i_sq_val, epsilon_0_val)
