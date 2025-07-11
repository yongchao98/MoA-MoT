import numpy as np
from scipy.integrate import trapz

def solve_density_profile():
    """
    Calculates the equilibrium mass density profile of a two-component, non-ideal gas
    mixture in a vertical cylinder under gravity using a self-consistent field method.
    """
    # 1. Define Constants and Parameters (in SI units)
    # Given parameters
    NA = 2.0e23  # number of particles of Gas A
    NB = 1.5e23 # number of particles of Gas B
    A = 0.1     # m^2, cross-sectional area
    H = 10.0    # m, height of container
    T = 500.0   # K, temperature
    g = 9.81    # m/s^2, gravitational acceleration

    # Molar masses (g/mol)
    M_A_mol = 28.0
    M_B_mol = 44.0

    # Van der Waals parameters (molar units: Pa * m^6 * mol^-2)
    a_AA_mol = 2.5
    a_BB_mol = 3.6
    a_AB_mol = 3.0

    # Physical constants
    N_avo = 6.02214076e23  # mol^-1, Avogadro constant
    k_B = 1.380649e-23     # J/K, Boltzmann constant

    # Convert parameters to per-particle SI units
    # Mass per particle (kg)
    m_A = M_A_mol / 1000.0 / N_avo
    m_B = M_B_mol / 1000.0 / N_avo

    # Interaction parameter 'a' (J*m^3)
    # Note: 1 Pa * m^6 = 1 J * m^3
    a_AA = a_AA_mol / N_avo**2
    a_BB = a_BB_mol / N_avo**2
    a_AB = a_AB_mol / N_avo**2

    # 2. Set up the Spatial Grid
    num_points = 201
    z = np.linspace(0, H, num_points)

    # 3. Iterative Solver
    # 3a. Initial guess (ideal gas barometric profiles)
    exp_factor_A_ideal = np.exp(-m_A * g * z / (k_B * T))
    exp_factor_B_ideal = np.exp(-m_B * g * z / (k_B * T))

    # Normalize to get the initial number density profiles
    integral_A_ideal = trapz(exp_factor_A_ideal, z) * A
    integral_B_ideal = trapz(exp_factor_B_ideal, z) * A
    n0_A_ideal = NA / integral_A_ideal
    n0_B_ideal = NB / integral_B_ideal

    n_A = n0_A_ideal * exp_factor_A_ideal
    n_B = n0_B_ideal * exp_factor_B_ideal

    # Iteration parameters
    num_iterations = 200
    mixing_factor = 0.1  # For stable convergence
    tolerance = 1e-12

    for i in range(num_iterations):
        n_A_old = n_A.copy()
        n_B_old = n_B.copy()

        # i. Calculate mean-field interaction potential
        U_int_A = -2 * (a_AA * n_A + a_AB * n_B)
        U_int_B = -2 * (a_BB * n_B + a_AB * n_A)

        # ii. Calculate total potential energy
        U_grav_A = m_A * g * z
        U_grav_B = m_B * g * z
        U_total_A = U_grav_A + U_int_A
        U_total_B = U_grav_B + U_int_B

        # iii. Calculate new unnormalized profiles via Boltzmann distribution
        # The potential at z=0 is subtracted to keep exponents manageable
        exp_factor_A = np.exp(-(U_total_A - U_total_A[0]) / (k_B * T))
        exp_factor_B = np.exp(-(U_total_B - U_total_B[0]) / (k_B * T))

        # iv. Normalize profiles to conserve total particle number
        integral_A = trapz(exp_factor_A, z) * A
        integral_B = trapz(exp_factor_B, z) * A
        
        # Avoid division by zero, though unlikely in this problem
        if integral_A == 0 or integral_B == 0:
            print("Error: Integral of profile is zero. Cannot normalize.")
            return

        n0_A = NA / integral_A
        n0_B = NB / integral_B
        
        n_A_new = n0_A * exp_factor_A
        n_B_new = n0_B * exp_factor_B

        # v. Mix new and old profiles for stability
        n_A = (1 - mixing_factor) * n_A_old + mixing_factor * n_A_new
        n_B = (1 - mixing_factor) * n_B_old + mixing_factor * n_B_new

        # Check for convergence
        change = np.linalg.norm(n_A - n_A_old) + np.linalg.norm(n_B - n_B_old)
        if change < tolerance:
            break

    # 4. Calculate Final Mass Density Profile
    rho_z = m_A * n_A + m_B * n_B

    # 5. Output the Results
    print("The final equation for the mass density profile is given by:")
    # Print the symbolic equation and values of the mass constants
    print(f"ρ(z) = m_A * n_A(z) + m_B * n_B(z)")
    print("\nWhere:")
    print(f"  m_A (mass of a Gas A particle) = {m_A:.4e} kg")
    print(f"  m_B (mass of a Gas B particle) = {m_B:.4e} kg")
    print("  n_A(z) and n_B(z) are the number density profiles determined numerically.")

    print("\nThe calculated mass density profile ρ(z) at various heights is:")
    print("-" * 40)
    print("Height z (m) | Mass Density ρ(z) (kg/m³)")
    print("-" * 40)

    # Print density at 11 points from z=0 to z=H
    for j in range(0, num_points, num_points // 10):
        print(f"{z[j]:>12.2f} | {rho_z[j]:>25.4f}")
        
    # Ensure the last point is always printed
    if (num_points - 1) % (num_points // 10) != 0:
        print(f"{z[-1]:>12.2f} | {rho_z[-1]:>25.4f}")
    print("-" * 40)

if __name__ == '__main__':
    solve_density_profile()