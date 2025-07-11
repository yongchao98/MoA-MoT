import numpy as np

def solve_density_profile():
    """
    Calculates the density profile of a non-ideal gas mixture in a gravitational field.
    """
    # --- Physical and System Constants ---
    N_A_total = 2e23  # Number of particles of Gas A
    N_B_total = 1.5e23 # Number of particles of Gas B
    M_A = 0.028  # Molar mass of Gas A (kg/mol)
    M_B = 0.044  # Molar mass of Gas B (kg/mol)
    A = 0.1  # Cross-sectional area of cylinder (m^2)
    H = 10.0  # Height of cylinder (m)
    T = 500.0  # Temperature (K)
    g = 9.81  # Gravitational acceleration (m/s^2)

    # Van der Waals parameters (molar)
    a_AA_molar = 2.5  # Pa * m^6 / mol^2
    b_AA_molar = 0.04  # m^3 / mol
    a_BB_molar = 3.6  # Pa * m^6 / mol^2
    b_BB_molar = 0.05  # m^3 / mol
    a_AB_molar = 3.0  # Pa * m^6 / mol^2

    # --- Fundamental Constants ---
    N_av = 6.02214076e23  # Avogadro's number (mol^-1)
    k_B = 1.380649e-23   # Boltzmann constant (J/K)
    
    # --- Convert Parameters to Per-Particle Basis ---
    m_A = M_A / N_av  # mass of a single particle A (kg)
    m_B = M_B / N_av  # mass of a single particle B (kg)

    a_AA = a_AA_molar / N_av**2
    a_BB = a_BB_molar / N_av**2
    a_AB = a_AB_molar / N_av**2
    b_AA = b_AA_molar / N_av
    b_BB = b_BB_molar / N_av
    
    # --- Numerical Simulation Parameters ---
    num_points = 201  # Number of points in the z-grid
    z = np.linspace(0, H, num_points) # z-grid from 0 to H
    
    max_iter = 2000    # Maximum number of iterations
    tolerance = 1e-8   # Convergence tolerance
    mix_alpha = 0.05   # Mixing factor for stability
    
    # --- Initialization ---
    V_total = A * H # Total Volume
    # Initial guess: uniform number density
    n_A = np.full(num_points, N_A_total / V_total)
    n_B = np.full(num_points, N_B_total / V_total)
    
    # --- Self-Consistent Field (SCF) Iteration ---
    for i in range(max_iter):
        n_A_old = n_A.copy()
        n_B_old = n_B.copy()
        
        # Calculate terms based on previous iteration's densities
        n_t = n_A + n_B  # Total number density
        n_b = n_A * b_AA + n_B * b_BB  # Excluded volume term
        
        # Avoid division by zero or log of non-positive by clamping n_b
        n_b = np.minimum(n_b, 0.9999999) 
        
        # Calculate the argument of the exponential for the Boltzmann-like factor
        # This comes from the chemical potential expression for a vdW mixture
        U_eff_A = (m_A * g * z) + (k_B * T * b_AA * n_t) / (1 - n_b) - 2 * (a_AA * n_A + a_AB * n_B)
        U_eff_B = (m_B * g * z) + (k_B * T * b_BB * n_t) / (1 - n_b) - 2 * (a_AB * n_A + a_BB * n_B)

        # Calculate unnormalized new profiles
        n_A_unnorm = (1 - n_b) * np.exp(-U_eff_A / (k_B * T))
        n_B_unnorm = (1 - n_b) * np.exp(-U_eff_B / (k_B * T))
        
        # Normalize profiles to conserve total particle number
        N_A_calc = A * np.trapz(n_A_unnorm, z)
        N_B_calc = A * np.trapz(n_B_unnorm, z)
        
        n_A_new = n_A_unnorm * (N_A_total / N_A_calc)
        n_B_new = n_B_unnorm * (N_B_total / N_B_calc)
        
        # Apply mixing to stabilize convergence
        n_A = mix_alpha * n_A_new + (1 - mix_alpha) * n_A_old
        n_B = mix_alpha * n_B_new + (1 - mix_alpha) * n_B_old
        
        # Check for convergence
        err_A = np.sqrt(np.mean((n_A - n_A_old)**2)) / np.mean(n_A)
        err_B = np.sqrt(np.mean((n_B - n_B_old)**2)) / np.mean(n_B)
        if err_A < tolerance and err_B < tolerance:
            # print(f"Converged after {i+1} iterations.")
            break

    # --- Final Calculation and Output ---
    # Calculate the total mass density profile
    rho_z = m_A * n_A + m_B * n_B
    
    print("The mass density profile rho(z) has been calculated.")
    print("-----------------------------------------------------")
    
    # Show the calculation for the density at the bottom of the container (z=0)
    print("Calculation for density at z = 0 meters:")
    print(f"rho(0) = m_A * n_A(0) + m_B * n_B(0)")
    print(f"rho(0) = ({m_A:.3e} kg) * ({n_A[0]:.3e} m^-3) + ({m_B:.3e} kg) * ({n_B[0]:.3e} m^-3)")
    print(f"rho(0) = {m_A * n_A[0]:.5f} kg/m^3 + {m_B * n_B[0]:.5f} kg/m^3 = {rho_z[0]:.5f} kg/m^3")
    print("-----------------------------------------------------")

    print("Mass density at selected heights:")
    # Print density at z=0, H/4, H/2, 3H/4, H
    height_indices = [0, num_points//4, num_points//2, 3*num_points//4, num_points-1]
    for idx in height_indices:
        print(f"rho(z = {z[idx]:>5.2f} m) = {rho_z[idx]:.5f} kg/m^3")
        
if __name__ == '__main__':
    solve_density_profile()
