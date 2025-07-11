import numpy as np

def solve_and_print_density_profile():
    """
    Calculates and prints the density profile of a non-ideal gas mixture
    in a gravitational field using a self-consistent mean-field model.
    """
    # --- Given Parameters ---
    N_A_particles = 2e23  # Number of particles of Gas A
    N_B_particles = 1.5e23 # Number of particles of Gas B
    M_A_molar = 28 / 1000  # Molar mass of Gas A in kg/mol
    M_B_molar = 44 / 1000  # Molar mass of Gas B in kg/mol
    A = 0.1  # Cross-sectional area in m^2
    H = 10.0 # Height of the container in m
    T = 500.0  # Temperature in K
    g = 9.81 # Gravitational acceleration in m/s^2
    a_AA_molar = 2.5  # van der Waals a for A-A in Pa*m^6*mol^-2
    a_BB_molar = 3.6  # van der Waals a for B-B in Pa*m^6*mol^-2
    a_AB_molar = 3.0  # van der Waals a for A-B in Pa*m^6*mol^-2

    # --- Physical Constants ---
    N_AVOGADRO = 6.02214076e23  # Avogadro's number in mol^-1
    k_B = 1.380649e-23  # Boltzmann constant in J/K

    # --- Derived Per-Particle Parameters (in SI units) ---
    m_A = M_A_molar / N_AVOGADRO  # Mass of a single particle of A in kg
    m_B = M_B_molar / N_AVOGADRO  # Mass of a single particle of B in kg
    
    # vdW 'a' parameters per particle pair. Units: J*m^3
    a_AA = a_AA_molar / N_AVOGADRO**2
    a_BB = a_BB_molar / N_AVOGADRO**2
    a_AB = a_AB_molar / N_AVOGADRO**2

    # --- Print the Governing Equations with Numerical Values ---
    print("The density profile ρ(z) is found by numerically solving the following self-consistent equations:")
    print("ρ_total_mass(z) = m_A * ρ_A(z) + m_B * ρ_B(z)\n")
    print("Where ρ_A(z) and ρ_B(z) are the number densities satisfying:")
    print("ρ_A(z) = C_A * exp(-(U_A(z)) / (k_B*T))")
    print("ρ_B(z) = C_B * exp(-(U_B(z)) / (k_B*T))\n")
    print("The total potential energy U(z) for each particle is:")
    print("U_A(z) = m_A*g*z - 2 * (a_AA*ρ_A(z) + a_AB*ρ_B(z))")
    print("U_B(z) = m_B*g*z - 2 * (a_BB*ρ_B(z) + a_AB*ρ_A(z))\n")
    print("Substituting the numerical values (in SI units):")
    print(f"k_B*T = {k_B*T:.4e} J")
    print(f"U_A(z) = ({m_A:.4e})*({g:.2f})*z - 2 * (({a_AA:.4e})*ρ_A(z) + ({a_AB:.4e})*ρ_B(z))")
    print(f"U_B(z) = ({m_B:.4e})*({g:.2f})*z - 2 * (({a_BB:.4e})*ρ_B(z) + ({a_AB:.4e})*ρ_A(z))")
    print("-" * 30)

    # --- Numerical Calculation Setup ---
    num_points = 200  # Number of points for discretization
    z_grid = np.linspace(0, H, num_points)
    dz = H / (num_points - 1)
    max_iter = 500  # Maximum number of iterations
    tolerance = 1e-12  # Convergence criterion
    mixing_factor = 0.1 # Mixing factor for stability

    # --- Initial Guess (Ideal Gas Barometric Formula) ---
    # Un-normalized profiles
    rho_A_unnormalized = np.exp(-m_A * g * z_grid / (k_B * T))
    rho_B_unnormalized = np.exp(-m_B * g * z_grid / (k_B * T))
    # Normalize to total particle count
    norm_A = np.sum(rho_A_unnormalized) * A * dz
    norm_B = np.sum(rho_B_unnormalized) * A * dz
    rho_A = (N_A_particles / norm_A) * rho_A_unnormalized
    rho_B = (N_B_particles / norm_B) * rho_B_unnormalized

    # --- Self-Consistent Iteration ---
    for i in range(max_iter):
        rho_A_old = rho_A.copy()
        rho_B_old = rho_B.copy()

        # Calculate interaction potential energy for each particle
        U_int_A = -2 * (a_AA * rho_A + a_AB * rho_B)
        U_int_B = -2 * (a_BB * rho_B + a_AB * rho_A)

        # Calculate total potential energy
        U_total_A = m_A * g * z_grid + U_int_A
        U_total_B = m_B * g * z_grid + U_int_B
        
        # To avoid overflow, subtract the minimum potential value before exponentiating
        U_total_A -= np.min(U_total_A)
        U_total_B -= np.min(U_total_B)

        # Update un-normalized profiles based on total potential
        rho_A_unnormalized_new = np.exp(-U_total_A / (k_B * T))
        rho_B_unnormalized_new = np.exp(-U_total_B / (k_B * T))

        # Normalize to get new profiles
        norm_A_new = np.sum(rho_A_unnormalized_new) * A * dz
        norm_B_new = np.sum(rho_B_unnormalized_new) * A * dz
        rho_A_new = (N_A_particles / norm_A_new) * rho_A_unnormalized_new
        rho_B_new = (N_B_particles / norm_B_new) * rho_B_unnormalized_new

        # Apply mixing for stable convergence
        rho_A = mixing_factor * rho_A_new + (1 - mixing_factor) * rho_A_old
        rho_B = mixing_factor * rho_B_new + (1 - mixing_factor) * rho_B_old

        # Check for convergence
        err = np.sum(np.abs(rho_A - rho_A_old)) + np.sum(np.abs(rho_B - rho_B_old))
        if err < tolerance:
            break

    # --- Calculate Final Mass Density Profile ---
    rho_mass_z = m_A * rho_A + m_B * rho_B

    # --- Output the final density profile at specific heights ---
    print("\nResulting Mass Density Profile ρ(z):")
    print(f"Density at z = 0.0 m (bottom): {rho_mass_z[0]:.5f} kg/m^3")
    print(f"Density at z = {H/2:.1f} m (middle): {rho_mass_z[num_points//2]:.5f} kg/m^3")
    print(f"Density at z = {H:.1f} m (top):   {rho_mass_z[-1]:.5f} kg/m^3")


if __name__ == "__main__":
    solve_and_print_density_profile()
