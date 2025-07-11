import numpy as np
from scipy.constants import g, Avogadro, k

def solve_and_print_density_profile():
    """
    Calculates and prints the density profile of a non-ideal gas mixture
    in a gravitational field using an iterative numerical method.
    """
    # 1. Parameters and Constants (in SI units)
    # System parameters
    A = 0.1  # Cross-sectional area (m^2)
    H = 10.0  # Height of container (m)
    T = 500.0  # Temperature (K)
    N_A = 2.0e23  # Number of particles of Gas A
    N_B = 1.5e23  # Number of particles of Gas B

    # Gas A parameters
    M_A = 28.0e-3  # Molar mass (kg/mol)
    a_AA = 2.5  # van der Waals 'a' (Pa * m^6 / mol^2)
    b_AA = 0.04e-3 # van der Waals 'b' (m^3 / mol) -> Corrected from prompt for consistency, typical values are in 10^-5 range, but using given value. Let's assume it's 0.04 m^3/mol, which is very large. Let's use a more physical 0.04 L/mol = 0.04e-3 m^3/mol. Let's stick to the prompt's values: 0.04 m^3/mol. This is likely a typo in the prompt as it's huge, but we follow it.
    b_AA = 0.04 # m^3/mol

    # Gas B parameters
    M_B = 44.0e-3  # Molar mass (kg/mol)
    a_BB = 3.6  # van der Waals 'a' (Pa * m^6 / mol^2)
    b_BB = 0.05 # van der Waals 'b' (m^3 / mol)

    # Interaction parameters
    a_AB = 3.0  # van der Waals 'a' for A-B interaction (Pa * m^6 / mol^2)

    # Convert molar parameters to per-particle parameters
    m_A = M_A / Avogadro  # mass of one particle of A (kg)
    m_B = M_B / Avogadro  # mass of one particle of B (kg)
    a_p_AA = a_AA / Avogadro**2
    a_p_BB = a_BB / Avogadro**2
    a_p_AB = a_AB / Avogadro**2
    b_p_AA = b_AA / Avogadro
    b_p_BB = b_BB / Avogadro
    
    k_B_T = k * T

    # Store parameters in a dictionary for easy access
    params = {
        'm_A': m_A, 'm_B': m_B, 'N_A': N_A, 'N_B': N_B,
        'a_p_AA': a_p_AA, 'a_p_BB': a_p_BB, 'a_p_AB': a_p_AB,
        'b_p_AA': b_p_AA, 'b_p_BB': b_p_BB, 'k_B_T': k_B_T, 'A': A
    }

    # 2. Numerical Settings
    z_points = 200
    z_grid = np.linspace(0, H, z_points)
    tolerance = 1e-8
    max_iter = 1000
    alpha = 0.05  # Mixing factor for stability

    # 3. Initial Guess (Ideal Gas Barometric Profile)
    # Calculate normalization constants for the ideal gas case
    beta_A = m_A * g / k_B_T
    beta_B = m_B * g / k_B_T
    integral_A = (1 - np.exp(-beta_A * H)) / beta_A
    integral_B = (1 - np.exp(-beta_B * H)) / beta_B
    n_A0 = N_A / (A * integral_A)
    n_B0 = N_B / (A * integral_B)
    
    n_A_current = n_A0 * np.exp(-beta_A * z_grid)
    n_B_current = n_B0 * np.exp(-beta_B * z_grid)

    # 4. Iterative Solver
    for i in range(max_iter):
        n_A_old = n_A_current.copy()
        n_B_old = n_B_current.copy()

        # Calculate excess chemical potential based on current density profiles
        vol_corr = 1.0 - (n_A_old * b_p_AA + n_B_old * b_p_BB)
        
        if np.any(vol_corr <= 0):
            print("Error: Unphysical density encountered (volume correction <= 0).")
            return

        log_term = -k_B_T * np.log(vol_corr)
        n_tot = n_A_old + n_B_old
        common_term = k_B_T * n_tot / vol_corr
        
        mu_ex_A = -2 * (n_A_old * a_p_AA + n_B_old * a_p_AB) + log_term + common_term * b_p_AA
        mu_ex_B = -2 * (n_A_old * a_p_AB + n_B_old * a_p_BB) + log_term + common_term * b_p_BB

        # Calculate new unnormalized profiles from equilibrium condition
        exponent_A = (-mu_ex_A - m_A * g * z_grid) / k_B_T
        exponent_B = (-mu_ex_B - m_B * g * z_grid) / k_B_T
        n_A_new_unnorm = np.exp(exponent_A)
        n_B_new_unnorm = np.exp(exponent_B)

        # Normalize to conserve total particle number
        N_A_calc = A * np.trapz(n_A_new_unnorm, z_grid)
        N_B_calc = A * np.trapz(n_B_new_unnorm, z_grid)
        
        C_A = N_A / N_A_calc
        C_B = N_B / N_B_calc

        n_A_new = C_A * n_A_new_unnorm
        n_B_new = C_B * n_B_new_unnorm

        # Mix new and old profiles for stable convergence
        n_A_current = alpha * n_A_new + (1 - alpha) * n_A_old
        n_B_current = alpha * n_B_new + (1 - alpha) * n_B_old

        # Check for convergence
        err_A = np.linalg.norm(n_A_current - n_A_old) / np.linalg.norm(n_A_old)
        err_B = np.linalg.norm(n_B_current - n_B_old) / np.linalg.norm(n_B_old)
        if err_A < tolerance and err_B < tolerance:
            print(f"Converged after {i+1} iterations.")
            break
    else:
        print("Warning: Maximum iterations reached without full convergence.")

    # 5. Final Calculation and Output
    rho_profile = n_A_current * m_A + n_B_current * m_B

    print("\n--- Density Profile rho(z) ---")
    print("The mass density profile rho(z) is determined numerically.")
    print("Mass density values (kg/m^3) at key heights:")
    
    idx_0 = 0
    idx_mid = len(z_grid) // 2
    idx_H = len(z_grid) - 1

    print(f"z = {z_grid[idx_0]:.1f} m: rho = {rho_profile[idx_0]:.4f} kg/m^3")
    print(f"z = {z_grid[idx_mid]:.1f} m: rho = {rho_profile[idx_mid]:.4f} kg/m^3")
    print(f"z = {z_grid[idx_H]:.1f} m: rho = {rho_profile[idx_H]:.4f} kg/m^3")

    print("\n--- Final Equation Details ---")
    print("The total mass density rho(z) is given by the equation:")
    print(f"rho(z) = n_A(z) * m_A + n_B(z) * m_B")
    print("\nwhere n_A(z) and n_B(z) are the number density profiles determined by this program.")
    print("The constant particle masses are:")
    print(f"m_A = M_A / N_av = {M_A:.4f} kg/mol / {Avogadro:.4e} mol^-1 = {m_A:.4e} kg")
    print(f"m_B = M_B / N_av = {M_B:.4f} kg/mol / {Avogadro:.4e} mol^-1 = {m_B:.4e} kg")
    
    print("\nUsing the calculated values at the middle of the container (z = 5.0 m) as an example:")
    print(f"n_A({z_grid[idx_mid]:.1f} m) = {n_A_current[idx_mid]:.4e} particles/m^3")
    print(f"n_B({z_grid[idx_mid]:.1f} m) = {n_B_current[idx_mid]:.4e} particles/m^3")
    print(f"rho({z_grid[idx_mid]:.1f}) = {n_A_current[idx_mid]:.4e} * {m_A:.4e} + {n_B_current[idx_mid]:.4e} * {m_B:.4e} = {rho_profile[idx_mid]:.4f} kg/m^3")

# Execute the main function
solve_and_print_density_profile()