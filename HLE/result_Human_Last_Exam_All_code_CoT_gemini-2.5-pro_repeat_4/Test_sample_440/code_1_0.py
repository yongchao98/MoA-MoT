import numpy as np

def solve_density_profile():
    """
    Calculates the density profile of a two-component non-ideal gas mixture
    in a gravitational field using a self-consistent field method.
    """
    # 1. Define Constants and Parameters
    # Physical constants
    N_AVOGADRO = 6.02214076e23  # mol^-1
    K_BOLTZMANN = 1.380649e-23  # J/K
    g = 9.81  # m/s^2

    # System parameters
    A = 0.1  # m^2
    H = 10.0  # m
    T = 500.0  # K
    N_A = 2.0e23  # particles
    N_B = 1.5e23  # particles

    # Gas properties (molar)
    M_A_molar = 28.0  # g/mol
    M_B_molar = 44.0  # g/mol
    a_AA_molar = 2.5  # Pa m^6 mol^-2
    a_BB_molar = 3.6  # Pa m^6 mol^-2
    a_AB_molar = 3.0  # Pa m^6 mol^-2

    # 2. Convert to Per-Particle SI Units
    # Mass per particle (kg)
    m_A = (M_A_molar / 1000) / N_AVOGADRO
    m_B = (M_B_molar / 1000) / N_AVOGADRO

    # Interaction parameter per particle pair (J m^3)
    # Note: 1 Pa = 1 J/m^3
    a_AA_p = a_AA_molar / N_AVOGADRO**2
    a_BB_p = a_BB_molar / N_AVOGADRO**2
    a_AB_p = a_AB_molar / N_AVOGADRO**2

    # 3. Setup Numerical Grid and Iteration Parameters
    num_points = 200
    z = np.linspace(0, H, num_points)
    dz = H / (num_points - 1)
    
    max_iter = 1000
    tolerance = 1e-9
    mix_factor = 0.1 # for stable convergence

    # 4. Initial Guess (Ideal Gas Barometric Distribution)
    # Calculate normalization constants for the ideal gas case
    integral_factor_A = (K_BOLTZMANN * T) / (m_A * g) * (1 - np.exp(-m_A * g * H / (K_BOLTZMANN * T)))
    n_A0_ideal = N_A / (A * integral_factor_A)
    integral_factor_B = (K_BOLTZMANN * T) / (m_B * g) * (1 - np.exp(-m_B * g * H / (K_BOLTZMANN * T)))
    n_B0_ideal = N_B / (A * integral_factor_B)
    
    n_A = n_A0_ideal * np.exp(-m_A * g * z / (K_BOLTZMANN * T))
    n_B = n_B0_ideal * np.exp(-m_B * g * z / (K_BOLTZMANN * T))

    # 5. Self-Consistent Iteration Loop
    for i in range(max_iter):
        n_A_old = n_A.copy()
        n_B_old = n_B.copy()

        # Calculate interaction potential energy for each species at each z
        U_int_A = -2 * (a_AA_p * n_A + a_AB_p * n_B)
        U_int_B = -2 * (a_BB_p * n_B + a_AB_p * n_A)
        
        # Calculate gravitational potential energy
        U_grav_A = m_A * g * z
        U_grav_B = m_B * g * z

        # Calculate the unnormalized new profiles (the exponential term)
        rhs_A = np.exp(-(U_grav_A + U_int_A) / (K_BOLTZMANN * T))
        rhs_B = np.exp(-(U_grav_B + U_int_B) / (K_BOLTZMANN * T))

        # Enforce particle number constraint to find normalization constants
        integral_A = np.sum(rhs_A) * dz * A
        integral_B = np.sum(rhs_B) * dz * A
        C_A = N_A / integral_A
        C_B = N_B / integral_B

        # Calculate new density profiles
        n_A_new = C_A * rhs_A
        n_B_new = C_B * rhs_B

        # Mix new and old profiles for stability
        n_A = mix_factor * n_A_new + (1 - mix_factor) * n_A_old
        n_B = mix_factor * n_B_new + (1 - mix_factor) * n_B_old
        
        # Check for convergence
        change = (np.sum(np.abs(n_A - n_A_old)) + np.sum(np.abs(n_B - n_B_old))) * dz
        if change < tolerance:
            # print(f"Converged after {i+1} iterations.")
            break

    # 6. Final Calculation and Output
    rho_z = m_A * n_A + m_B * n_B

    print("--- Model and Parameters ---")
    print("The mass density profile rho(z) is calculated based on a self-consistent mean-field model.")
    print("The number densities n_A(z) and n_B(z) solve the coupled equations:")
    print("n_A(z) = C_A * exp(-(m_A*g*z - 2*(a'_AA*n_A(z) + a'_AB*n_B(z))) / (k_B*T))")
    print("n_B(z) = C_B * exp(-(m_B*g*z - 2*(a'_BB*n_B(z) + a'_AB*n_A(z))) / (k_B*T))")
    print("\nwhere rho(z) = m_A*n_A(z) + m_B*n_B(z).")
    print("\nKey parameter values (SI units):")
    print(f"m_A (mass of particle A): {m_A:.4e} kg")
    print(f"m_B (mass of particle B): {m_B:.4e} kg")
    print(f"a'_AA (interaction parameter): {a_AA_p:.4e} J*m^3")
    print(f"a'_BB (interaction parameter): {a_BB_p:.4e} J*m^3")
    print(f"a'_AB (interaction parameter): {a_AB_p:.4e} J*m^3")
    print(f"k_B (Boltzmann constant): {K_BOLTZMANN:.4e} J/K")
    print(f"T (Temperature): {T:.1f} K")
    print(f"g (Gravitational acceleration): {g:.2f} m/s^2")

    print("\n--- Resulting Density Profile ---")
    print(f"Mass density at the bottom (z=0 m): {rho_z[0]:.6f} kg/m^3")
    print(f"Mass density at the middle (z={H/2:.1f} m): {rho_z[num_points//2]:.6f} kg/m^3")
    print(f"Mass density at the top (z={H:.1f} m): {rho_z[-1]:.6f} kg/m^3")

    # Final answer as requested by the format
    print(f"\n<<<{rho_z[0]:.6f}>>>")

if __name__ == '__main__':
    solve_density_profile()