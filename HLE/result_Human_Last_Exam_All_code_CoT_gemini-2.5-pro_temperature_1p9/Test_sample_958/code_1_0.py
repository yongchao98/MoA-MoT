import numpy as np

def calculate_stored_energy_coeffs(n, theta_deg):
    """
    Calculates the numerical coefficients for the stored energy expressions
    from option D, using given n and theta.
    Constants like omega/c, epsilon_0, and E_x0_i are set to 1 for simplicity.
    """
    # Convert angle to radians
    theta_rad = np.deg2rad(theta_deg)
    
    # Check for TIR condition
    theta_c_rad = np.arcsin(1/n)
    if theta_rad <= theta_c_rad:
        print(f"Error: Total Internal Reflection requires theta > theta_c.")
        print(f"For n={n}, critical angle theta_c is {np.rad2deg(theta_c_rad):.2f} degrees.")
        print(f"Provided theta = {theta_deg} degrees is not greater than the critical angle.")
        return None, None

    # Intermediate calculations
    n2 = n**2
    s2 = np.sin(theta_rad)**2
    
    # Common terms in the denominator
    # Let w_c = omega/c = 1
    # Let eps0 = 1
    # Let E_x0_i_sq = 1
    
    term_sqrt = np.sqrt(n2 * s2 - 1)
    term_denom1 = n2 - 1
    term_denom2 = (n2 + 1) * s2 - 1
    
    denominator = term_denom1 * term_denom2 * term_sqrt
    
    # Numerator for Electric field energy
    num_E = n2 * (2 * n2 * s2 - 1)
    
    # Numerator for Magnetic field energy
    num_H = n2 * (n2 * s2 - 1)
    
    # Calculate final coefficients
    coeff_E = num_E / denominator
    coeff_H = num_H / denominator
    
    return coeff_E, coeff_H, num_E, num_H, term_denom1, term_denom2, term_sqrt

# --- Parameters ---
# Refractive index of the core material
n_core = 1.5
# Angle of incidence in degrees
theta_incidence_deg = 60

print("This script calculates the numerical coefficients for the time-averaged stored energy per unit area based on the formulas in Option D.\n")
print("We use the following example values:")
print(f"Refractive index n = {n_core}")
print(f"Angle of incidence theta = {theta_incidence_deg} degrees")
print("For simplicity, we assume (omega/c) = 1, epsilon_0 = 1, and |E_{x0}^i|^2 = 1.\n")


# Perform the calculation
coeffs = calculate_stored_energy_coeffs(n_core, theta_incidence_deg)

if coeffs:
    coeff_E, coeff_H, num_E, num_H, d1, d2, d_sqrt = coeffs
    
    print("Based on Option D, the formulas are:")
    print("Energy in E field = (n^2 * (2*n^2*sin^2(theta) - 1)) / ((omega/c) * (n^2-1) * ((n^2+1)*sin^2(theta)-1) * sqrt(n^2*sin^2(theta)-1)) * epsilon_0 * |E_{x0}^i|^2")
    print("Energy in H field = (n^2 * (n^2*sin^2(theta) - 1)) / ((omega/c) * (n^2-1) * ((n^2+1)*sin^2(theta)-1) * sqrt(n^2*sin^2(theta)-1)) * epsilon_0 * |E_{x0}^i|^2")
    print("-" * 20)

    # Print the detailed breakdown of the calculation as requested
    print("Calculation breakdown for the dimensionless coefficient:")
    print("\nFor the Electric Field Energy:")
    print(f"  Numerator = n^2 * (2*n^2*sin^2(theta) - 1) = {n_core**2:.4f} * (2*{n_core**2:.4f}*{np.sin(np.deg2rad(theta_incidence_deg))**2:.4f} - 1) = {num_E:.4f}")
    print(f"  Denominator = (n^2-1) * ((n^2+1)*sin^2(theta)-1) * sqrt(n^2*sin^2(theta)-1)")
    print(f"              = {d1:.4f} * {d2:.4f} * {d_sqrt:.4f} = {d1*d2*d_sqrt:.4f}")
    print(f"Resulting Coefficient for E-field energy: {coeff_E:.4f}")
    
    print("\nFor the Magnetic Field Energy:")
    print(f"  Numerator = n^2 * (n^2*sin^2(theta) - 1) = {n_core**2:.4f} * ({n_core**2:.4f}*{np.sin(np.deg2rad(theta_incidence_deg))**2:.4f} - 1) = {num_H:.4f}")
    print(f"  Denominator is the same: {d1*d2*d_sqrt:.4f}")
    print(f"Resulting Coefficient for H-field energy: {coeff_H:.4f}")
