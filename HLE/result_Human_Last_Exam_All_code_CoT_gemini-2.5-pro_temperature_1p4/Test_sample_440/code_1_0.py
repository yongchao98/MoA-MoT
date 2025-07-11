import numpy as np

def solve_density_profile():
    """
    This function calculates the density profile of a two-component, non-ideal gas
    mixture in a gravitational field using a first-order perturbation method.
    """
    # --- 1. Constants and Parameters (in SI units) ---
    N_A = 2.0e23
    N_B = 1.5e23
    M_A_mol = 28e-3  # kg/mol
    M_B_mol = 44e-3  # kg/mol
    A = 0.1
    H = 10.0
    T = 500.0
    g = 9.81
    a_AA_mol = 2.5
    a_BB_mol = 3.6
    a_AB_mol = 3.0
    N_avo = 6.02214076e23
    k_B = 1.380649e-23

    # Derived per-particle parameters
    m_A = M_A_mol / N_avo
    m_B = M_B_mol / N_avo
    a_AA = a_AA_mol / N_avo**2
    a_BB = a_BB_mol / N_avo**2
    a_AB = a_AB_mol / N_avo**2

    # Numerical grid
    num_points = 1001
    z = np.linspace(0, H, num_points)
    dz = H / (num_points - 1)
    kBT = k_B * T

    # --- 2. Zeroth-Order (Ideal Gas) Solution ---
    integral_exp_A = (kBT / (m_A * g)) * (1 - np.exp(-m_A * g * H / kBT))
    integral_exp_B = (kBT / (m_B * g)) * (1 - np.exp(-m_B * g * H / kBT))
    C_A0 = N_A / (A * integral_exp_A)
    C_B0 = N_B / (A * integral_exp_B)
    n_A0_z = C_A0 * np.exp(-m_A * g * z / kBT)
    n_B0_z = C_B0 * np.exp(-m_B * g * z / kBT)

    # --- 3. First-Order Correction (Interaction Potential) ---
    U_int_A_z = -2 * (a_AA * n_A0_z + a_AB * n_B0_z)
    U_int_B_z = -2 * (a_BB * n_B0_z + a_AB * n_A0_z)

    # --- 4. Unnormalized First-Order Density Profiles ---
    potential_A_z = (m_A * g * z + U_int_A_z) / kBT
    potential_B_z = (m_B * g * z + U_int_B_z) / kBT
    n_A_unnorm_z = np.exp(-potential_A_z)
    n_B_unnorm_z = np.exp(-potential_B_z)

    # --- 5. Renormalization ---
    integral_unnorm_A = np.sum(n_A_unnorm_z) * dz
    integral_unnorm_B = np.sum(n_B_unnorm_z) * dz
    Norm_A = N_A / (A * integral_unnorm_A)
    Norm_B = N_B / (A * integral_unnorm_B)
    n_A_z = Norm_A * n_A_unnorm_z
    n_B_z = Norm_B * n_B_unnorm_z

    # --- 6. Final Mass Density Profile ---
    rho_z = m_A * n_A_z + m_B * n_B_z

    # --- 7. Output Results ---
    idx_H_half = np.argmin(np.abs(z - H / 2.0))
    n_A_mid = n_A_z[idx_H_half]
    n_B_mid = n_B_z[idx_H_half]
    rho_A_mid = m_A * n_A_mid
    rho_B_mid = m_B * n_B_mid
    rho_total_mid = rho_z[idx_H_half]

    print("--- Detailed Calculation for Density at Midpoint z = 5.0 m ---")
    print("The final equation for density is: rho(z) = m_A * n_A(z) + m_B * n_B(z)\n")
    print(f"The numbers for the equation at z = 5.0 m are:")
    print(f"rho_A = m_A * n_A = {m_A:.3e} kg * {n_A_mid:.3e} m^-3 = {rho_A_mid:.5f} kg/m^3")
    print(f"rho_B = m_B * n_B = {m_B:.3e} kg * {n_B_mid:.3e} m^-3 = {rho_B_mid:.5f} kg/m^3")
    print(f"rho_total = rho_A + rho_B = {rho_A_mid:.5f} kg/m^3 + {rho_B_mid:.5f} kg/m^3 = {rho_total_mid:.5f} kg/m^3")

solve_density_profile()