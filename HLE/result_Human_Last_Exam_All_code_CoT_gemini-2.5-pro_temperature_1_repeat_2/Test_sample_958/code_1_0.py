import numpy as np

def calculate_stored_energy(n, theta_deg, omega_c=1.0, epsilon0=1.0, Ex0_i_sq=1.0):
    """
    Calculates the time-averaged stored energy per unit area in the electric and
    magnetic fields of an evanescent wave for p-polarization.

    Args:
        n (float): Refractive index of the core.
        theta_deg (float): Incident angle in degrees.
        omega_c (float): Angular frequency divided by the speed of light (k0).
        epsilon0 (float): Permittivity of free space.
        Ex0_i_sq (float): Squared magnitude of the x-component of the incident electric field.
    """
    # Convert angle to radians
    theta = np.deg2rad(theta_deg)

    # Check for TIR condition
    theta_c_rad = np.arcsin(1/n)
    if theta <= theta_c_rad:
        print("Warning: Angle is not greater than the critical angle. Evanescent wave is not generated under these conditions.")
        print(f"Critical angle: {np.rad2deg(theta_c_rad):.2f} degrees")
        return None, None

    sin2_theta = np.sin(theta)**2
    n2 = n**2

    # Common denominator part of the expressions
    term1 = n2 * sin2_theta - 1
    if term1 < 0:
        # This should not happen if theta > theta_c
        print("Error: n*sin(theta) < 1. Cannot calculate square root of a negative number.")
        return None, None
        
    sqrt_term = np.sqrt(term1)
    denom_part = (n2 - 1) * ((n2 + 1) * sin2_theta - 1) * sqrt_term
    common_factor = 1 / (2 * omega_c * denom_part)

    # Calculate energy in the Electric field
    e_field_num = n2 * (2 * n2 * sin2_theta - 1)
    energy_E = e_field_num * common_factor * epsilon0 * Ex0_i_sq

    # Calculate energy in the Magnetic field
    h_field_num = n2 * (n2 * sin2_theta - 1)
    energy_H = h_field_num * common_factor * epsilon0 * Ex0_i_sq
    
    # Printing the results
    print("For the given parameters:")
    print(f"  Refractive index n = {n}")
    print(f"  Incident angle theta = {theta_deg} degrees")
    print(f"  (w/c) = {omega_c}, epsilon_0 = {epsilon0}, |E_x0^i|^2 = {Ex0_i_sq}\n")
    
    print("The final equations are:")
    print("Energy in E field = (n^2 * (2*n^2*sin^2(theta) - 1)) / (2 * (w/c) * (n^2 - 1) * ((n^2 + 1)*sin^2(theta) - 1) * sqrt(n^2*sin^2(theta) - 1)) * epsilon_0 * |E_x0^i|^2")
    print("Energy in H field = (n^2 * (n^2*sin^2(theta) - 1)) / (2 * (w/c) * (n^2 - 1) * ((n^2 + 1)*sin^2(theta) - 1) * sqrt(n^2*sin^2(theta) - 1)) * epsilon_0 * |E_x0^i|^2\n")

    print("Calculated values:")
    print(f"Energy in E field = {energy_E}")
    print(f"Energy in H field = {energy_H}")

    return energy_E, energy_H

# Example usage with typical values
# Glass core (n=1.5), air cladding (n=1)
# Incident angle must be greater than the critical angle arcsin(1/1.5) = 41.8 degrees.
n_core = 1.5
angle_incident = 60.0
calculate_stored_energy(n=n_core, theta_deg=angle_incident)