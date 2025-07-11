import numpy as np

def solve_density_profile():
    """
    Calculates the density profile of a two-component non-ideal gas mixture
    in a gravitational field using a self-consistent field (SCF) method.
    """
    # --- 1. Define Constants and Parameters (in SI units) ---
    # System parameters
    A = 0.1  # Cross-sectional area (m^2)
    H = 10.0  # Height of container (m)
    T = 500.0  # Temperature (K)
    g = 9.81  # Gravitational acceleration (m/s^2)
    
    # Gas A properties
    N_A = 2.0e23  # Number of particles of Gas A
    M_A = 0.028   # Molar mass of Gas A (kg/mol)
    a_AA = 2.5    # vdW parameter a_AA (Pa * m^6 * mol^-2)
    
    # Gas B properties
    N_B = 1.5e23  # Number of particles of Gas B
    M_B = 0.044   # Molar mass of Gas B (kg/mol)
    a_BB = 3.6    # vdW parameter a_BB (Pa * m^6 * mol^-2)

    # Interaction parameter
    a_AB = 3.0    # vdW interaction parameter a_AB (Pa * m^6 * mol^-2)

    # Physical constants
    N_avogadro = 6.02214076e23  # Avogadro's number (mol^-1)
    R = 8.314462618           # Ideal gas constant (J * mol^-1 * K^-1)
    k_B = R / N_avogadro      # Boltzmann constant (J/K)

    # --- 2. Pre-calculations ---
    # Particle masses
    m_A = M_A / N_avogadro
    m_B = M_B / N_avogadro

    # Convert vdW 'a' parameters from molar to per-particle basis
    # a' = a / N_avogadro^2, units: J*m^3
    a_prime_AA = a_AA / N_avogadro**2
    a_prime_BB = a_BB / N_avogadro**2
    a_prime_AB = a_AB / N_avogadro**2

    # --- 3. Setup Numerical Grid and Iteration Parameters ---
    Nz = 200  # Number of grid points for height
    z = np.linspace(0, H, Nz)
    dz = H / (Nz - 1)
    
    num_iterations = 200 # Number of SCF iterations
    mixing_factor = 0.1  # Mixing factor for stable convergence

    # --- 4. Initialize Density Profiles (using ideal gas barometric formula) ---
    # This provides a good starting point for the SCF iteration.
    # n_i(0) = N_i / integral(A * exp(-m_i*g*z/(k_B*T)) dz)
    term_A = m_A * g * H / (k_B * T)
    integral_factor_A = (A * k_B * T / (m_A * g)) * (1 - np.exp(-term_A))
    n_A0_ideal = N_A / integral_factor_A
    n_A_z = n_A0_ideal * np.exp(-m_A * g * z / (k_B * T))

    term_B = m_B * g * H / (k_B * T)
    integral_factor_B = (A * k_B * T / (m_B * g)) * (1 - np.exp(-term_B))
    n_B0_ideal = N_B / integral_factor_B
    n_B_z = n_B0_ideal * np.exp(-m_B * g * z / (k_B * T))

    # --- 5. Self-Consistent Field (SCF) Iteration ---
    for i in range(num_iterations):
        n_A_old = n_A_z.copy()
        n_B_old = n_B_z.copy()

        # Calculate mean-field interaction potential energy for each particle type
        # U_int_A = -2 * (a'_AA*n_A + a'_AB*n_B)
        U_int_A = -2 * (a_prime_AA * n_A_old + a_prime_AB * n_B_old)
        U_int_B = -2 * (a_prime_AB * n_A_old + a_prime_BB * n_B_old)
        
        # Calculate new unnormalized profiles based on total potential
        # (gravitational + interaction)
        n_A_unnorm = np.exp(-(m_A * g * z + U_int_A) / (k_B * T))
        n_B_unnorm = np.exp(-(m_B * g * z + U_int_B) / (k_B * T))
        
        # Normalize profiles to conserve total particle number
        integral_A = np.trapz(n_A_unnorm, z) * A
        integral_B = np.trapz(n_B_unnorm, z) * A
        
        C_A = N_A / integral_A
        C_B = N_B / integral_B
        
        n_A_new = C_A * n_A_unnorm
        n_B_new = C_B * n_B_unnorm
        
        # Apply mixing to stabilize convergence
        n_A_z = mixing_factor * n_A_new + (1 - mixing_factor) * n_A_old
        n_B_z = mixing_factor * n_B_new + (1 - mixing_factor) * n_B_old

    # --- 6. Final Calculations and Output ---
    # Total mass density profile
    rho_z = m_A * n_A_z + m_B * n_B_z

    print("The density profile ρ(z) is determined by the following set of self-consistent equations:")
    print("-" * 80)
    
    # Output the equations with all numerical values plugged in
    print(f"ρ(z) = ({m_A:.3e}) * n_A(z) + ({m_B:.3e}) * n_B(z)  [kg/m^3]\n")
    print("where n_A(z) and n_B(z) are solved from the coupled equations:\n")
    
    # Equation for n_A(z)
    term1_A = f"({m_A * g:.3e} * z"
    term2_A = f" - 2 * (({a_prime_AA:.3e})*n_A(z) + ({a_prime_AB:.3e})*n_B(z))"
    denom = f") / {k_B * T:.3e}"
    print(f"n_A(z) = {C_A:.3e} * exp( -({term1_A}{term2_A}{denom}) )")

    # Equation for n_B(z)
    term1_B = f"({m_B * g:.3e} * z"
    term2_B = f" - 2 * (({a_prime_AB:.3e})*n_A(z) + ({a_prime_BB:.3e})*n_B(z))"
    print(f"n_B(z) = {C_B:.3e} * exp( -({term1_B}{term2_B}{denom}) )\n")
    
    print("-" * 80)
    print("The numerical solution to these equations yields the following density values:")
    
    # Find indices for z = 0, H/2, and H
    idx_bottom = 0
    idx_middle = np.argmin(np.abs(z - H / 2.0))
    idx_top = Nz - 1
    
    print(f"Density at the bottom (z = 0 m):   ρ(0)   = {rho_z[idx_bottom]:.4f} kg/m^3")
    print(f"Density at the middle (z = {H/2} m):  ρ(H/2) = {rho_z[idx_middle]:.4f} kg/m^3")
    print(f"Density at the top (z = {H} m):    ρ(H)   = {rho_z[idx_top]:.4f} kg/m^3")

if __name__ == '__main__':
    solve_density_profile()