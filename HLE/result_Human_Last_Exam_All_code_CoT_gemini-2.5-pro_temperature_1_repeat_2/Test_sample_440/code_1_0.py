import numpy as np

def solve_density_profile():
    """
    Calculates the density profile of a two-component non-ideal gas mixture
    in a gravitational field using a self-consistent field iteration.
    """
    # --- Physical and System Constants ---
    # Given parameters
    A = 0.1  # Cross-sectional area (m^2)
    H = 10.0  # Height of container (m)
    T = 500.0  # Temperature (K)
    g = 9.81  # Gravitational acceleration (m/s^2)
    N_A = 2.0e23  # Number of particles of Gas A
    N_B = 1.5e23  # Number of particles of Gas B
    M_A_gmol = 28.0  # Molar mass of Gas A (g/mol)
    M_B_gmol = 44.0  # Molar mass of Gas B (g/mol)
    a_AA_molar = 2.5  # van der Waals a_AA (Pa m^6 mol^-2 or J m^3 mol^-2)
    a_BB_molar = 3.6  # van der Waals a_BB (Pa m^6 mol^-2 or J m^3 mol^-2)
    a_AB_molar = 3.0  # van der Waals a_AB (Pa m^6 mol^-2 or J m^3 mol^-2)

    # Fundamental constants
    N_avo = 6.02214076e23  # Avogadro's number (mol^-1)
    k_B = 1.380649e-23  # Boltzmann constant (J/K)

    # --- Derived Per-Particle Parameters ---
    # Convert molar mass from g/mol to kg/particle
    m_A = (M_A_gmol / 1000) / N_avo  # Mass of one particle of Gas A (kg)
    m_B = (M_B_gmol / 1000) / N_avo  # Mass of one particle of Gas B (kg)

    # Convert molar vdW 'a' to per-particle 'a_prime' (J m^3)
    a_prime_AA = a_AA_molar / N_avo**2
    a_prime_BB = a_BB_molar / N_avo**2
    a_prime_AB = a_AB_molar / N_avo**2

    # --- Numerical Setup ---
    num_points = 200  # Number of points for discretization of height
    z = np.linspace(0, H, num_points)
    dz = H / (num_points - 1)
    iterations = 100 # Number of self-consistent iterations
    mixing_alpha = 0.1 # Mixing factor for convergence stability

    # --- Initial Guess (Ideal Gas Barometric Formula) ---
    # Un-normalized profiles
    n_A_unnorm = np.exp(-m_A * g * z / (k_B * T))
    n_B_unnorm = np.exp(-m_B * g * z / (k_B * T))

    # Normalization constants C_A and C_B
    integral_A = A * np.trapz(n_A_unnorm, z)
    integral_B = A * np.trapz(n_B_unnorm, z)
    C_A = N_A / integral_A
    C_B = N_B / integral_B

    # Initial number density profiles (particles/m^3)
    n_A = C_A * n_A_unnorm
    n_B = C_B * n_B_unnorm
    
    # --- Self-Consistent Field Iteration ---
    for i in range(iterations):
        # Calculate mean-field interaction potential from vdW 'a' term
        U_int_A = -2 * (a_prime_AA * n_A + a_prime_AB * n_B)
        U_int_B = -2 * (a_prime_AB * n_A + a_prime_BB * n_B)

        # Calculate new un-normalized profiles including interaction potential
        n_A_unnorm_new = np.exp(-(m_A * g * z + U_int_A) / (k_B * T))
        n_B_unnorm_new = np.exp(-(m_B * g * z + U_int_B) / (k_B * T))

        # Normalize the new profiles to conserve particle number
        integral_A_new = A * np.trapz(n_A_unnorm_new, z)
        integral_B_new = A * np.trapz(n_B_unnorm_new, z)
        C_A_new = N_A / integral_A_new
        C_B_new = N_B / integral_B_new
        
        n_A_new = C_A_new * n_A_unnorm_new
        n_B_new = C_B_new * n_B_unnorm_new

        # Mix old and new profiles for stable convergence
        n_A = (1 - mixing_alpha) * n_A + mixing_alpha * n_A_new
        n_B = (1 - mixing_alpha) * n_B + mixing_alpha * n_B_new

    # --- Final Results ---
    # The final converged normalization constants
    C_A_final = N_A / (A * np.trapz(np.exp(-(m_A * g * z - 2 * (a_prime_AA * n_A + a_prime_AB * n_B)) / (k_B * T)), z))
    C_B_final = N_B / (A * np.trapz(np.exp(-(m_B * g * z - 2 * (a_prime_AB * n_A + a_prime_BB * n_B)) / (k_B * T)), z))

    # Calculate final mass density profile (kg/m^3)
    rho_z = n_A * m_A + n_B * m_B

    # --- Print Output ---
    print("The density profile rho(z) is determined by the following self-consistent equations:")
    print("-" * 70)
    print(f"rho(z) = ({m_A:.3e} kg) * n_A(z) + ({m_B:.3e} kg) * n_B(z)\n")
    
    print("where n_A(z) and n_B(z) are the number densities (particles/m^3) satisfying:")
    print(f"n_A(z) = C_A * exp( -[({m_A*g:.3e} J/m)*z - 2*(({a_prime_AA:.3e} J m^3)*n_A(z) + ({a_prime_AB:.3e} J m^3)*n_B(z))] / ({k_B*T:.3e} J) )")
    print(f"n_B(z) = C_B * exp( -[({m_B*g:.3e} J/m)*z - 2*(({a_prime_AB:.3e} J m^3)*n_A(z) + ({a_prime_BB:.3e} J m^3)*n_B(z))] / ({k_B*T:.3e} J) )\n")
    
    print("The numerically determined normalization constants are:")
    print(f"C_A ≈ {C_A_final:.3e} m^-3")
    print(f"C_B ≈ {C_B_final:.3e} m^-3\n")
    
    print("-" * 70)
    print("The resulting mass density at the bottom of the container is:")
    print(f"rho(z=0) = {rho_z[0]:.4f} kg/m^3")
    print("-" * 70)
    
    return rho_z[0]

# Execute the function and capture the final answer
final_answer = solve_density_profile()
print(f"<<<{final_answer:.4f}>>>")