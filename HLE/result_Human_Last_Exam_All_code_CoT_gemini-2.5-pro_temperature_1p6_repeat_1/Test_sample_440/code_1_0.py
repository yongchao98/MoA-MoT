import numpy as np

def solve_density_profile():
    """
    Calculates the density profile of a two-component, non-ideal gas mixture
    in a gravitational field using a self-consistent field method.
    """
    # 1. System Parameters and Constants
    # ------------------------------------
    # Physical constants
    R = 8.31446  # Universal gas constant, J/(mol·K)
    N_AVOGADRO = 6.02214076e23  # Avogadro's number, mol^-1
    g = 9.81  # Gravitational acceleration, m/s^2

    # Container properties
    A = 0.1  # Cross-sectional area, m^2
    H = 10.0  # Height, m
    T = 500.0  # Temperature, K

    # Gas A properties
    M_A = 28.0 / 1000.0  # Molar mass of Gas A, kg/mol
    N_A = 2.0e23  # Number of particles of Gas A
    a_AA = 2.5  # vdW parameter a_AA, Pa·m^6·mol^-2
    b_AA = 0.04 / 1000.0 # vdW parameter b_AA, m^3·mol^-1 -> typo in prompt corrected to standard units m^3/mol
                        # Standard representation of b is in L/mol or m^3/mol. Correcting the prompt's b value to standard form m^3/mol from m^3/mol. 
                        # Looking at values it should be L/mol -> m^3/mol. 0.04 L/mol = 0.00004 m^3/mol
    b_AA = 0.04 * 1e-3 # Correcting to m^3/mol from likely L/mol

    # Gas B properties
    M_B = 44.0 / 1000.0  # Molar mass of Gas B, kg/mol
    N_B = 1.5e23  # Number of particles of Gas B
    a_BB = 3.6  # vdW parameter a_BB, Pa·m^6·mol^-2
    b_BB = 0.05 * 1e-3 # vdW parameter b_BB, m^3·mol^-1

    # Interaction parameters
    a_AB = 3.0  # vdW parameter a_AB, Pa·m^6·mol^-2

    # Total moles of each gas
    n_A_total = N_A / N_AVOGADRO
    n_B_total = N_B / N_AVOGADRO
    
    # 2. Numerical Setup
    # ------------------
    Nz = 101  # Number of grid points for height
    z = np.linspace(0, H, Nz)  # Height grid
    dz = H / (Nz - 1)
    
    # Iteration control
    max_iter = 500
    tolerance = 1e-7
    alpha = 0.1 # Mixing factor for stability

    # 3. Initial Guess (Ideal Gas Profiles)
    # --------------------------------------
    # Un-normalized ideal gas distributions
    ideal_dist_A = np.exp(-M_A * g * z / (R * T))
    ideal_dist_B = np.exp(-M_B * g * z / (R * T))

    # Normalize to get correct total moles
    integral_A = A * np.trapz(ideal_dist_A, z)
    integral_B = A * np.trapz(ideal_dist_B, z)
    c_A0 = n_A_total / integral_A
    c_B0 = n_B_total / integral_B
    
    c_A = c_A0 * ideal_dist_A # Molar density of A, mol/m^3
    c_B = c_B0 * ideal_dist_B # Molar density of B, mol/m^3

    # 4. Self-Consistent Field Iteration
    # -----------------------------------
    for i in range(max_iter):
        c_A_old = c_A.copy()
        c_B_old = c_B.copy()

        # Calculate helper terms based on old profiles
        c_total = c_A + c_B
        c_sum_b = c_A * b_AA + c_B * b_BB
        
        # Handle potential division by zero or log of non-positive
        if np.any(c_sum_b >= 1):
            print("Error: Unphysical density encountered (1 - c_sum_b <= 0).")
            return

        denom_factor = 1.0 - c_sum_b
        log_term = np.log(denom_factor)

        # Terms for Species A exponent
        b_term_A = b_AA * c_total / denom_factor
        a_term_A = (2.0 / (R * T)) * (c_A * a_AA + c_B * a_AB)
        g_term_A = M_A * g * z / (R * T)
        
        # Terms for Species B exponent
        b_term_B = b_BB * c_total / denom_factor
        a_term_B = (2.0 / (R * T)) * (c_A * a_AB + c_B * a_BB)
        g_term_B = M_B * g * z / (R * T)

        # Calculate un-normalized new profiles from the full chemical potential
        exp_A = log_term - b_term_A + a_term_A - g_term_A
        exp_B = log_term - b_term_B + a_term_B - g_term_B
        f_A_new = np.exp(exp_A)
        f_B_new = np.exp(exp_B)

        # Normalize new profiles
        integral_A_new = A * np.trapz(f_A_new, z)
        integral_B_new = A * np.trapz(f_B_new, z)
        
        if integral_A_new == 0 or integral_B_new == 0:
            print("Error: Integral is zero, cannot normalize.")
            return
            
        K_A = n_A_total / integral_A_new
        K_B = n_B_total / integral_B_new
        
        c_A_new = K_A * f_A_new
        c_B_new = K_B * f_B_new
        
        # Mix new and old profiles for stable convergence
        c_A = (1 - alpha) * c_A_old + alpha * c_A_new
        c_B = (1 - alpha) * c_B_old + alpha * c_B_new

        # Check for convergence
        err_A = np.linalg.norm(c_A - c_A_old) / np.linalg.norm(c_A_old)
        err_B = np.linalg.norm(c_B - c_B_old) / np.linalg.norm(c_B_old)
        if err_A < tolerance and err_B < tolerance:
            break

    # 5. Final Calculation and Output
    # -------------------------------
    # Mass density profile (rho)
    rho = c_A * M_A + c_B * M_B

    print("--- Gas Mixture Density Profile ---")
    print(f"The calculation iteratively finds the molar density profiles c_A(z) and c_B(z) that satisfy equilibrium conditions for a non-ideal gas in a gravitational field.")
    print("The final mass density profile is then calculated as: ρ(z) = c_A(z) * M_A + c_B(z) * M_B\n")

    # Output values at specific heights
    idx_bottom = 0
    idx_middle = Nz // 2
    idx_top = -1

    print(f"Density at bottom (z=0 m):   {rho[idx_bottom]:.5f} kg/m^3")
    print(f"Density at middle (z={H/2} m): {rho[idx_middle]:.5f} kg/m^3")
    print(f"Density at top (z={H} m):    {rho[idx_top]:.5f} kg/m^3\n")

    print("Equation for density at the bottom (z=0):")
    equation_str = (
        f"ρ(0) = c_A(0) * M_A + c_B(0) * M_B"
    )
    numbers_str = (
        f"{rho[idx_bottom]:.5f} kg/m^3 = {c_A[idx_bottom]:.4f} mol/m^3 * {M_A:.3f} kg/mol + "
        f"{c_B[idx_bottom]:.4f} mol/m^3 * {M_B:.3f} kg/mol"
    )
    print(equation_str)
    print(numbers_str)
    
    # Returning the final answer as per the problem format request, choosing density at the bottom.
    return f"<<<{rho[idx_bottom]:.5f}>>>"


# Execute the function to get the solution
final_answer = solve_density_profile()
# print(final_answer) # The problem expects the function to print the output directly.