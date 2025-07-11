import numpy as np

def solve_gas_density_profile():
    """
    Calculates the density profile of a non-ideal two-gas mixture
    in a vertical cylinder under gravity.
    """
    # --- Parameters ---
    # Physical constants
    R = 8.314  # J/(mol*K)
    N_AVOGADRO = 6.022e23  # particles/mol
    g = 9.81  # m/s^2

    # Container properties
    A = 0.1  # m^2
    H = 10.0  # m
    T = 500.0  # K

    # Gas A properties
    N_A = 2e23  # number of particles
    M_A_g_mol = 28.0  # g/mol
    M_A = M_A_g_mol / 1000.0  # kg/mol

    # Gas B properties
    N_B = 1.5e23  # number of particles
    M_B_g_mol = 44.0  # g/mol
    M_B = M_B_g_mol / 1000.0  # kg/mol

    # Van der Waals interaction parameters
    a_AA = 2.5  # Pa * m^6 / mol^2
    a_BB = 3.6  # Pa * m^6 / mol^2
    a_AB = 3.0  # Pa * m^6 / mol^2

    # --- Derived quantities ---
    n_A_total = N_A / N_AVOGADRO  # total moles of A
    n_B_total = N_B / N_AVOGADRO  # total moles of B

    # --- Numerical parameters ---
    num_points = 201
    z = np.linspace(0, H, num_points)
    max_iterations = 200
    tolerance = 1e-7

    # --- Initialization: Start with Ideal Gas Profile ---
    MG_A = M_A * g
    MG_B = M_B * g
    RT = R * T
    
    integral_factor_A = (RT / MG_A) * (1 - np.exp(-MG_A * H / RT)) if MG_A > 0 else H
    integral_factor_B = (RT / MG_B) * (1 - np.exp(-MG_B * H / RT)) if MG_B > 0 else H

    c_A0_ideal = n_A_total / (A * integral_factor_A)
    c_B0_ideal = n_B_total / (A * integral_factor_B)
    
    # Initial guess for molar density profiles (mol/m^3)
    c_A = c_A0_ideal * np.exp(-MG_A * z / RT)
    c_B = c_B0_ideal * np.exp(-MG_B * z / RT)

    # --- Iterative Solution ---
    for i in range(max_iterations):
        c_A_old = c_A.copy()
        c_B_old = c_B.copy()

        # Calculate the mean-field interaction potential energy per mole
        V_int_A = -2 * (c_A * a_AA + c_B * a_AB)
        V_int_B = -2 * (c_A * a_AB + c_B * a_BB)

        # Total effective potential per mole
        V_eff_A = MG_A * z + V_int_A
        V_eff_B = MG_B * z + V_int_B

        # Calculate un-normalized profiles, shifting potential to prevent overflow
        f_A = np.exp(-(V_eff_A - np.min(V_eff_A)) / RT)
        f_B = np.exp(-(V_eff_B - np.min(V_eff_B)) / RT)
        
        # Re-normalize to conserve total moles
        integral_A = np.trapz(f_A, z)
        integral_B = np.trapz(f_B, z)
        
        C_prime_A = n_A_total / (A * integral_A)
        C_prime_B = n_B_total / (A * integral_B)
        
        c_A = C_prime_A * f_A
        c_B = C_prime_B * f_B
        
        # Check for convergence
        error = np.max(np.abs(c_A - c_A_old))
        if error < tolerance:
            break

    # --- Final Calculation ---
    # Calculate the total mass density profile (kg/m^3)
    rho_z = c_A * M_A + c_B * M_B

    # --- Output Results ---
    print("Density Profile of the Gas Mixture ρ(z)\n")
    print("The mass density ρ(z) is calculated from the converged molar densities c_A(z) and c_B(z) and their respective molar masses M_A and M_B.")
    print(f"M_A = {M_A:.4f} kg/mol, M_B = {M_B:.4f} kg/mol\n")

    # Find indices for z = 0, H/2, and H
    idx_bottom = 0
    idx_middle = np.argmin(np.abs(z - H / 2.0))
    idx_top = num_points - 1
    indices_to_print = [idx_bottom, idx_middle, idx_top]

    for i in indices_to_print:
        z_val = z[i]
        c_A_val = c_A[i]
        c_B_val = c_B[i]
        rho_val = rho_z[i]
        print(f"At height z = {z_val:.2f} m:")
        print(f"  ρ(z) = c_A(z) * M_A + c_B(z) * M_B")
        print(f"       = {c_A_val:.4f} mol/m^3 * {M_A:.4f} kg/mol + {c_B_val:.4f} mol/m^3 * {M_B:.4f} kg/mol")
        print(f"       = {rho_val:.4f} kg/m^3\n")
        
    return rho_z[0]

# Execute the function and get the final answer
final_answer = solve_gas_density_profile()
print(f"<<<{final_answer}>>>")