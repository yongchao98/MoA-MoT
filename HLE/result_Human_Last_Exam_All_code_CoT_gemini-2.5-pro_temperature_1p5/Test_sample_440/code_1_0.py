import numpy as np
from scipy.integrate import quad

def solve_density_profile():
    """
    Calculates the density profile of a two-component, non-ideal gas mixture
    in a gravitational field using a first-order correction to the barometric formula.
    """
    # 1. Define Constants and Parameters
    # System parameters
    A = 0.1  # Cross-sectional area [m^2]
    H = 10   # Height of the container [m]
    T = 500  # Temperature [K]
    g = 9.81 # Gravitational acceleration [m/s^2]
    N_A_particles = 2e23
    N_B_particles = 1.5e23

    # Gas A parameters
    M_A_gmol = 28 # g/mol
    a_AA = 2.5    # Pa * m^6 * mol^-2
    b_AA = 0.04   # m^3 * mol^-1 (not used in this model)

    # Gas B parameters
    M_B_gmol = 44 # g/mol
    a_BB = 3.6    # Pa * m^6 * mol^-2
    b_BB = 0.05   # m^3 * mol^-1 (not used in this model)
    
    # Interaction parameter
    a_AB = 3.0 # Pa * m^6 * mol^-2

    # Physical constants
    N_avo = 6.02214076e23 # Avogadro's number [mol^-1]
    R = 8.314462618     # Ideal gas constant [J/(mol*K)]

    # Convert molar masses to SI units [kg/mol]
    M_A = M_A_gmol / 1000
    M_B = M_B_gmol / 1000

    # Calculate total moles of each gas
    n_A_tot = N_A_particles / N_avo
    n_B_tot = N_B_particles / N_avo
    
    # --- Step 1: Ideal Gas Profile Calculation (Zeroth-Order) ---
    
    # Integrand for the normalization integral of the barometric formula
    def integrand_ideal(z, M, g_const, R_const, T_const):
        return np.exp(-M * g_const * z / (R_const * T_const))

    # Perform the normalization integrals
    integral_norm_A_ideal, _ = quad(integrand_ideal, 0, H, args=(M_A, g, R, T))
    integral_norm_B_ideal, _ = quad(integrand_ideal, 0, H, args=(M_B, g, R, T))

    # Calculate the ideal gas molar densities at z=0
    c_A0_ideal = n_A_tot / (A * integral_norm_A_ideal)
    c_B0_ideal = n_B_tot / (A * integral_norm_B_ideal)

    # Define the ideal molar density profile functions
    def c_A_ideal(z):
        return c_A0_ideal * np.exp(-M_A * g * z / (R * T))

    def c_B_ideal(z):
        return c_B0_ideal * np.exp(-M_B * g * z / (R * T))

    # --- Step 2: Mean-Field Potential and Corrected Profile ---
    
    # Define interaction potential energy based on the ideal gas densities
    def U_int_A(z):
        return -2 * (a_AA * c_A_ideal(z) + a_AB * c_B_ideal(z))

    def U_int_B(z):
        return -2 * (a_AB * c_A_ideal(z) + a_BB * c_B_ideal(z))

    # Define integrands for the new normalization integrals (proportional to final densities)
    def integrand_corrected_A(z):
        U_total_A = M_A * g * z + U_int_A(z)
        return np.exp(-U_total_A / (R * T))

    def integrand_corrected_B(z):
        U_total_B = M_B * g * z + U_int_B(z)
        return np.exp(-U_total_B / (R * T))

    # Calculate the new normalization integrals
    integral_norm_A_corr, _ = quad(integrand_corrected_A, 0, H)
    integral_norm_B_corr, _ = quad(integrand_corrected_B, 0, H)

    # Calculate the normalization constants for the corrected densities
    C_A = n_A_tot / (A * integral_norm_A_corr)
    C_B = n_B_tot / (A * integral_norm_B_corr)

    # --- Step 3: Final Density Profile Functions ---
    
    def c_A_final(z):
        return C_A * integrand_corrected_A(z)

    def c_B_final(z):
        return C_B * integrand_corrected_B(z)

    def rho_final(z):
        return c_A_final(z) * M_A + c_B_final(z) * M_B

    # --- Step 4: Output the Results ---
    
    print("The final mass density profile is given by the function rho(z).")
    print("rho(z) = c_A_final(z) * M_A + c_B_final(z) * M_B")
    print("\n--- Key Parameters and Constants ---")
    print(f"M_A = {M_A:.4f} kg/mol")
    print(f"M_B = {M_B:.4f} kg/mol")
    print("c_A_final(z) = C_A * exp(-(M_A*g*z + U_int_A(z)) / (R*T))")
    print("c_B_final(z) = C_B * exp(-(M_B*g*z + U_int_B(z)) / (R*T))")
    print(f"Calculated Normalization Constant C_A = {C_A:.4f} mol/m^3")
    print(f"Calculated Normalization Constant C_B = {C_B:.4f} mol/m^3")
    print(f"(Note: U_int(z) is calculated using zeroth-order ideal gas densities)")

    print("\n--- Density Profile rho(z) [kg/m^3] ---")
    z_points = [0, H / 2, H]
    for z_val in z_points:
        density = rho_final(z_val)
        print(f"rho(z = {z_val:2.1f} m) = {density:.5f}")
        
    # As a representative value, calculate the average density
    total_mass = n_A_tot * M_A + n_B_tot * M_B
    total_volume = A * H
    average_density = total_mass / total_volume
    print(f"\nFor reference, the average density is: {average_density:.5f} kg/m^3")


if __name__ == '__main__':
    solve_density_profile()
