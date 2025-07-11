import numpy as np
from scipy.integrate import quad

def solve_density_profile():
    """
    Calculates and describes the density profile of a non-ideal gas mixture
    in a gravitational field.
    """
    # 1. Parameters and Constants
    N_A_particles = 2e23
    N_B_particles = 1.5e23
    M_A = 0.028  # kg/mol
    M_B = 0.044  # kg/mol
    A = 0.1      # m^2
    H = 10       # m
    T = 500      # K
    g = 9.81     # m/s^2
    a_AA = 2.5   # Pa m^6 mol^-2, which is J m^3 mol^-2
    a_BB = 3.6   # J m^3 mol^-2
    a_AB = 3.0   # J m^3 mol^-2
    R = 8.314    # J/(mol K)
    N_Av = 6.022e23 # mol^-1

    # Convert particle numbers to moles
    n_A_total = N_A_particles / N_Av
    n_B_total = N_B_particles / N_Av

    # 2. Ideal Gas Reference Profile Calculation
    RT = R * T
    # Gravitational exponents
    alpha_A = M_A * g / RT
    alpha_B = M_B * g / RT

    # Integrate to find normalization constants for the ideal case
    # The integral of exp(-alpha*z) from 0 to H is (1 - exp(-alpha*H))/alpha
    integral_A_ideal = (1 - np.exp(-alpha_A * H)) / alpha_A
    integral_B_ideal = (1 - np.exp(-alpha_B * H)) / alpha_B

    # Molar density at z=0 for the ideal case
    c_A0_ideal = n_A_total / (A * integral_A_ideal)
    c_B0_ideal = n_B_total / (A * integral_B_ideal)

    # Define ideal gas profile functions
    def c_A_ideal(z):
        return c_A0_ideal * np.exp(-alpha_A * z)
    def c_B_ideal(z):
        return c_B0_ideal * np.exp(-alpha_B * z)

    # 3. Corrected Density Profile Calculation
    # Define integrands for the corrected profiles (unnormalized)
    def integrand_A_corr(z):
        # Interaction potential energy per mole for Gas A at height z
        U_int_A = -2 * (a_AA * c_A_ideal(z) + a_AB * c_B_ideal(z))
        # Total potential energy per mole
        U_total_A = M_A * g * z + U_int_A
        return np.exp(-U_total_A / RT)

    def integrand_B_corr(z):
        # Interaction potential energy per mole for Gas B at height z
        U_int_B = -2 * (a_BB * c_B_ideal(z) + a_AB * c_A_ideal(z))
        # Total potential energy per mole
        U_total_B = M_B * g * z + U_int_B
        return np.exp(-U_total_B / RT)

    # Numerically integrate to find the new normalization constants
    integral_A_corr, _ = quad(integrand_A_corr, 0, H)
    integral_B_corr, _ = quad(integrand_B_corr, 0, H)

    # Normalization constants for the corrected profiles
    C_A_corr = n_A_total / (A * integral_A_corr)
    C_B_corr = n_B_total / (A * integral_B_corr)

    # Define final, corrected molar density profile functions
    def c_A_corr(z):
        return C_A_corr * integrand_A_corr(z)

    def c_B_corr(z):
        return C_B_corr * integrand_B_corr(z)

    # 4. Total Mass Density Profile
    def total_mass_density(z):
        """Calculates the total mass density rho(z) in kg/m^3."""
        return c_A_corr(z) * M_A + c_B_corr(z) * M_B

    # 5. Output the Results
    print("DETERMINATION OF THE GAS MIXTURE DENSITY PROFILE ρ(z)")
    print("====================================================\n")
    print("This model determines the density profile using a Boltzmann distribution corrected")
    print("for mean-field molecular attractions (van der Waals 'a' parameter). The excluded")
    print("volume effect ('b' parameter) is neglected. The calculation uses a perturbative")
    print("approach where the interaction potential is estimated from an ideal gas density profile.\n")
    print(f"The final mass density profile is given by ρ(z) = {M_A} * c_A(z) + {M_B} * c_B(z).\n")

    print("--- Component A Profile: c_A(z) [mol/m^3] ---")
    print(f"c_A(z) = C_A_corr * exp(-(M_A*g*z + U_int_A(z)) / (R*T))")
    print("The components of this equation are:")
    print(f"  - C_A_corr (Normalization Constant) = {C_A_corr:.5f} mol/m^3")
    print("  - Interaction Potential U_int_A(z) = -2 * (a_AA*c_A_ideal(z) + a_AB*c_B_ideal(z))")
    print(f"    - a_AA = {a_AA} J*m^3/mol^2")
    print(f"    - a_AB = {a_AB} J*m^3/mol^2")
    print("  - Ideal Gas Reference c_A_ideal(z) = c_A0_ideal * exp(-M_A*g*z / (R*T))")
    print(f"    - c_A0_ideal = {c_A0_ideal:.5f} mol/m^3")
    print("\n--- Component B Profile: c_B(z) [mol/m^3] ---")
    print(f"c_B(z) = C_B_corr * exp(-(M_B*g*z + U_int_B(z)) / (R*T))")
    print("The components of this equation are:")
    print(f"  - C_B_corr (Normalization Constant) = {C_B_corr:.5f} mol/m^3")
    print("  - Interaction Potential U_int_B(z) = -2 * (a_BB*c_B_ideal(z) + a_AB*c_A_ideal(z))")
    print(f"    - a_BB = {a_BB} J*m^3/mol^2")
    print(f"    - a_AB = {a_AB} J*m^3/mol^2")
    print("  - Ideal Gas Reference c_B_ideal(z) = c_B0_ideal * exp(-M_B*g*z / (R*T))")
    print(f"    - c_B0_ideal = {c_B0_ideal:.5f} mol/m^3")

    print("\n--- Total Mass Density ρ(z) at Specific Heights ---")
    z_points = [0, H/2, H]
    for z in z_points:
        rho_val = total_mass_density(z)
        print(f"ρ(z={z:<4.1f} m) = {rho_val:.5f} kg/m^3")
        
# Execute the function to print the solution
solve_density_profile()