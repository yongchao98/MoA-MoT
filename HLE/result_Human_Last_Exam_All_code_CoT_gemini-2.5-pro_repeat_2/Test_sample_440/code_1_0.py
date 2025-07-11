import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root

def solve_gas_profile():
    """
    Solves for the density profile of a two-component non-ideal gas mixture
    in a gravitational field using the van der Waals equation of state.
    """
    # --- Parameters (all in SI units) ---
    # Gas A
    M_A = 28e-3  # Molar mass of Gas A (kg/mol)
    N_A = 2e23   # Number of particles of Gas A
    # Gas B
    M_B = 44e-3  # Molar mass of Gas B (kg/mol)
    N_B = 1.5e23 # Number of particles of Gas B
    # van der Waals parameters
    a_AA = 2.5   # Pa * m^6 * mol^-2
    b_AA = 0.04  # m^3 * mol^-1 (Note: This value is unusually large, but used as given)
    a_BB = 3.6   # Pa * m^6 * mol^-2
    b_BB = 0.05  # m^3 * mol^-1 (Note: This value is unusually large, but used as given)
    a_AB = 3.0   # Pa * m^6 * mol^-2
    # Container
    A = 0.1      # Cross-sectional area (m^2)
    H = 10.0     # Height (m)
    # Physical constants
    T = 500.0    # Temperature (K)
    g = 9.81     # Gravitational acceleration (m/s^2)
    R = 8.314    # Ideal gas constant (J * mol^-1 * K^-1)
    N_AVOGADRO = 6.022e23 # Avogadro's number (mol^-1)

    # Target total moles for each gas
    moles_A_target = N_A / N_AVOGADRO
    moles_B_target = N_B / N_AVOGADRO

    def get_derivatives(z, n_vec):
        """
        Calculates the derivatives dn_A/dz and dn_B/dz at a given height z.
        This function defines the system of ordinary differential equations.
        n_vec: [n_A, n_B], vector of molar densities (mol/m^3).
        """
        n_A, n_B = n_vec

        if n_A <= 1e-9 or n_B <= 1e-9: # Prevent non-physical or unstable values
            return [0, 0]

        # Intermediate terms for the Jacobian of scaled chemical potentials
        n_t = n_A + n_B
        beta = n_A * b_AA + n_B * b_BB # Total excluded volume fraction term
        D = 1.0 - beta

        if D <= 1e-9: # Non-physical state (excluded volume exceeds total volume)
            return [np.inf, np.inf] # Signal error to the ODE solver

        # Jacobian elements J_ij = d(mu_i/RT)/dn_j
        J11 = (1/n_A) + (b_AA/D) + b_AA * (1 - beta + n_t * b_AA) / (D**2) - (2 * a_AA) / (R * T)
        J12 = (b_BB/D) + b_AA * (1 - beta + n_t * b_BB) / (D**2) - (2 * a_AB) / (R * T)
        J21 = (b_AA/D) + b_BB * (1 - beta + n_t * b_AA) / (D**2) - (2 * a_AB) / (R * T)
        J22 = (1/n_B) + (b_BB/D) + b_BB * (1 - beta + n_t * b_BB) / (D**2) - (2 * a_BB) / (R * T)

        J = np.array([[J11, J12], [J21, J22]])
        rhs = -np.array([M_A * g, M_B * g]) / (R * T)

        try:
            dndt = np.linalg.solve(J, rhs)
        except np.linalg.LinAlgError:
            return [np.inf, np.inf] # Signal error if Jacobian is singular

        return dndt

    def objective_function(n0_vec):
        """
        Objective function for the root finder (shooting method).
        Calculates the error between calculated and target total moles for a given n(0).
        """
        if np.any(n0_vec <= 0): # Initial guess must be physical
            return [1e6, 1e6]
            
        sol = solve_ivp(get_derivatives, [0, H], n0_vec, method='RK45', dense_output=True)

        if sol.status != 0: # Check for integration failure
            return [1e6, 1e6]
        
        z_nodes = np.linspace(0, H, 201)
        n_profiles = sol.sol(z_nodes)
        
        if np.any(n_profiles < 0): # Check for non-physical results
             return [1e6, 1e6]

        moles_A_calc = np.trapz(n_profiles[0], z_nodes) * A
        moles_B_calc = np.trapz(n_profiles[1], z_nodes) * A

        return [moles_A_calc - moles_A_target, moles_B_calc - moles_B_target]

    print("Starting the numerical calculation for the density profile...")
    
    # Initial guess for densities at z=0 (use average density)
    V_total = A * H
    initial_guess = [moles_A_target / V_total, moles_B_target / V_total]

    # Find the correct initial densities n(0)
    solution = root(objective_function, initial_guess, method='hybr', tol=1e-8)

    if not solution.success:
        print("\nError: The numerical solver failed to find the correct initial conditions.")
        print(f"Solver message: {solution.message}")
        return

    n0_correct = solution.x
    print("Solver successfully found the initial conditions at z=0.")

    # Solve the ODE one last time with the correct n(0) to get the final profiles
    z_eval = np.linspace(0, H, 11)
    final_sol = solve_ivp(get_derivatives, [0, H], n0_correct, method='RK45', t_eval=z_eval)

    n_A_z = final_sol.y[0]
    n_B_z = final_sol.y[1]
    
    # Calculate the total mass density profile rho(z)
    rho_z = n_A_z * M_A + n_B_z * M_B

    print("\n" + "="*50)
    print("<<< Gas Mixture Density Profile ρ(z) >>>")
    print("="*50)
    print("\nThe density profile ρ(z) is determined by solving the hydrostatic\nequilibrium conditions for a non-ideal van der Waals gas mixture.")
    print("\nBelow are the calculated mass densities at various heights:")
    print("-" * 42)
    print(f"{'Height z (m)':^18} | {'Density ρ(z) (kg/m^3)':^22}")
    print("-" * 42)
    for z, rho in zip(z_eval, rho_z):
        print(f"{z:^18.2f} | {rho:^22.4f}")
    print("-" * 42)

    # Print the equation for the first point as requested
    print("\nThe density at any height z is given by the equation:")
    print("  ρ(z) = n_A(z) * M_A + n_B(z) * M_B")
    print("\nExample calculation for z = 0.00 m:")
    print(f"  ρ(0) = n_A(0) * M_A + n_B(0) * M_B")
    print(f"  ρ(0) = {n_A_z[0]:.4f} [mol/m^3] * {M_A} [kg/mol] + {n_B_z[0]:.4f} [mol/m^3] * {M_B} [kg/mol]")
    print(f"  ρ(0) = {(n_A_z[0] * M_A):.4f} [kg/m^3] + {(n_B_z[0] * M_B):.4f} [kg/m^3]")
    print(f"  ρ(0) = {rho_z[0]:.4f} [kg/m^3]")
    print("="*50)


if __name__ == '__main__':
    solve_gas_profile()