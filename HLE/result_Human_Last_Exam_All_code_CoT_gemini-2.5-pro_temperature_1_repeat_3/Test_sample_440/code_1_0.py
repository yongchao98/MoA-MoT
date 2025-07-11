import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import trapz

def solve_density_profile():
    """
    Calculates and prints the density profile of a non-ideal gas mixture
    in a gravitational field.
    """
    # 1. Define constants and parameters in SI units
    # Given parameters
    N_A_total = 2e23      # particles
    N_B_total = 1.5e23    # particles
    M_A_molar = 0.028     # kg/mol
    M_B_molar = 0.044     # kg/mol
    A = 0.1               # m^2
    H = 10                # m
    T = 500               # K
    g = 9.81              # m/s^2
    a_AA_molar = 2.5      # Pa m^6 mol^-2
    b_AA_molar = 0.04     # m^3 mol^-1
    a_BB_molar = 3.6      # Pa m^6 mol^-2
    b_BB_molar = 0.05     # m^3 mol^-1
    a_AB_molar = 3.0      # Pa m^6 mol^-2

    # Physical constants
    N_av = 6.02214076e23  # mol^-1
    k_B = 1.380649e-23    # J/K

    # Convert parameters to per-particle units for calculations
    m_A = M_A_molar / N_av  # kg/particle
    m_B = M_B_molar / N_av  # kg/particle
    a_AA = a_AA_molar / N_av**2
    a_BB = a_BB_molar / N_av**2
    a_AB = a_AB_molar / N_av**2
    b_A = b_AA_molar / N_av   # m^3/particle
    b_B = b_BB_molar / N_av   # m^3/particle
    
    # 2. Define helper function for chemical potential
    def get_chemical_potential(n_A, n_B):
        """Calculates the intrinsic chemical potential for each species."""
        vol_frac = n_A * b_A + n_B * b_B
        if n_A <= 0 or n_B <= 0 or vol_frac >= 1:
            return np.inf, np.inf # Unphysical state

        log_term_n_A = k_B * T * np.log(n_A)
        log_term_n_B = k_B * T * np.log(n_B)
        
        common_log_term = -k_B * T * np.log(1 - vol_frac)
        b_term_denom = 1 / (1 - vol_frac)

        mu_A_prime = common_log_term + k_B * T * n_A * b_A * b_term_denom - 2 * (n_A * a_AA + n_B * a_AB)
        mu_B_prime = common_log_term + k_B * T * n_B * b_B * b_term_denom - 2 * (n_A * a_AB + n_B * a_BB)

        return log_term_n_A + mu_A_prime, log_term_n_B + mu_B_prime

    # Discretize height for numerical integration
    z_grid = np.linspace(0, H, 101)

    # 3. Define the main function for the root-finding algorithm
    memo_profiles = {}
    def get_profiles_and_error(n0_guess):
        """
        For a guess of densities at z=0, calculates the full density profiles
        and returns the error in total particle numbers.
        """
        n_A0, n_B0 = n0_guess
        
        # Use memoization to avoid re-computing for the same guess
        if tuple(n0_guess) in memo_profiles:
            return memo_profiles[tuple(n0_guess)]

        # Calculate reference total chemical potential at z=0
        mu_A_0, mu_B_0 = get_chemical_potential(n_A0, n_B0)
        if np.isinf(mu_A_0): # Unphysical guess
            return (1e30, 1e30), (None, None)

        # Arrays to store density profiles
        n_A_profile = np.zeros_like(z_grid)
        n_B_profile = np.zeros_like(z_grid)
        n_A_profile[0], n_B_profile[0] = n_A0, n_B0

        # Define equations to solve for n(z) at a given height z
        def equations_at_z(n_guess, z):
            mu_A, mu_B = get_chemical_potential(n_guess[0], n_guess[1])
            if np.isinf(mu_A):
                return [1e30, 1e30]
            # Error is the difference from the constant total potential
            err_A = mu_A + m_A * g * z - mu_A_0
            err_B = mu_B + m_B * g * z - mu_B_0
            return [err_A, err_B]

        # Solve for n(z) at each grid point
        current_n_guess = [n_A0, n_B0]
        for i, z in enumerate(z_grid[1:], 1):
            solution, _, ier, _ = fsolve(equations_at_z, current_n_guess, args=(z,), full_output=True)
            if ier != 1: # Solver failed
                return (1e30, 1e30), (None, None)
            n_A_profile[i], n_B_profile[i] = solution
            current_n_guess = solution
            
        # Integrate profiles to get total particle numbers
        N_A_calc = trapz(n_A_profile, z_grid) * A
        N_B_calc = trapz(n_B_profile, z_grid) * A
        
        error = ((N_A_calc - N_A_total), (N_B_calc - N_B_total))
        profiles = (n_A_profile, n_B_profile)
        memo_profiles[tuple(n0_guess)] = error, profiles
        return error, profiles

    def error_function_for_solver(n0_guess):
        """Wrapper function for the outer solver."""
        error, _ = get_profiles_and_error(n0_guess)
        return error

    # 4. Solve for the base densities n(0)
    # Initial guess: average density, slightly increased to account for gravity
    V_total = A * H
    n_A_avg = N_A_total / V_total
    n_B_avg = N_B_total / V_total
    initial_guess_n0 = [1.2 * n_A_avg, 1.2 * n_B_avg]

    n0_solution, _, ier, msg = fsolve(error_function_for_solver, initial_guess_n0, full_output=True)
    if ier != 1:
        print("Warning: Solver for base densities may not have converged.")
        print("Message:", msg)
        
    # 5. Calculate final profiles with the correct base densities
    _, (n_A_profile, n_B_profile) = get_profiles_and_error(n0_solution)
    rho_profile = n_A_profile * m_A + n_B_profile * m_B
    mu_A_const, mu_B_const = get_chemical_potential(n0_solution[0], n0_solution[1])

    # 6. Print the results
    print("--- Model and Parameters ---")
    print("The density profile rho(z) is determined by solving the equilibrium conditions:")
    print("  μ_A(n_A(z), n_B(z)) + m_A*g*z = C_A")
    print("  μ_B(n_A(z), n_B(z)) + m_B*g*z = C_B")
    print("where μ_i is the van der Waals chemical potential and n_i(z) is the number density profile.")
    print("\nKey parameters (SI units):")
    print(f"  Particle mass A (m_A): {m_A:.3e} kg")
    print(f"  Particle mass B (m_B): {m_B:.3e} kg")
    print(f"  Gravitational acceleration (g): {g} m/s^2")
    print(f"  Temperature (T): {T} K")
    print(f"  Boltzmann constant (k_B): {k_B:.3e} J/K")
    print(f"  vdW params (a_AA, a_AB, a_BB): {a_AA:.3e}, {a_AB:.3e}, {a_BB:.3e}")
    print(f"  vdW params (b_A, b_B): {b_A:.3e}, {b_B:.3e}")
    print("\nNumerically determined constants:")
    print(f"  Total chemical potential C_A: {mu_A_const:.3e} J")
    print(f"  Total chemical potential C_B: {mu_B_const:.3e} J")
    print("\n--- Density Profile rho(z) ---")
    print(" z (m) | rho(z) (kg/m^3)")
    print("----------------------------")
    # Print profile at 11 points from z=0 to z=H
    for i in range(0, len(z_grid), len(z_grid)//10):
        z_val = z_grid[i]
        rho_val = rho_profile[i]
        print(f" {z_val:5.2f} | {rho_val:15.5f}")

if __name__ == '__main__':
    solve_density_profile()
