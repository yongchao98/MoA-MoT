import numpy as np
from scipy.integrate import solve_ivp, trapz
from scipy.optimize import root

def solve_density_profile():
    """
    Calculates the density profile of a two-component non-ideal gas mixture
    in a gravitational field using a numerical shooting method.
    """
    # 1. Define Constants and Parameters (all in SI units)
    A = 0.1         # Cross-sectional area (m^2)
    H = 10.0        # Height of container (m)
    T = 500.0       # Temperature (K)
    g = 9.81        # Gravitational acceleration (m/s^2)
    N_A_particles = 2e23  # Number of particles of Gas A
    N_B_particles = 1.5e23 # Number of particles of Gas B

    # Molar masses
    M_A = 28e-3     # kg/mol
    M_B = 44e-3     # kg/mol

    # Van der Waals parameters
    a_AA = 2.5      # Pa m^6 mol^-2
    a_BB = 3.6      # Pa m^6 mol^-2
    a_AB = 3.0      # Pa m^6 mol^-2
    b_AA = 0.04     # m^3 mol^-1
    b_BB = 0.05     # m^3 mol^-1

    # Physical constants
    R = 8.31446     # J/(mol K)
    N_avo = 6.022e23 # mol^-1

    # Total moles of each gas
    n_tot_A = N_A_particles / N_avo
    n_tot_B = N_B_particles / N_avo

    # 2. Define the system of ODEs: d(n_i)/dz
    def get_derivatives(z, n):
        """Calculates dn/dz for the given molar densities n = [nA, nB]."""
        nA, nB = n[0], n[1]

        # Ensure physical validity (densities must be positive)
        if nA <= 0 or nB <= 0:
            return [np.inf, np.inf]
        
        # Excluded volume term from van der Waals equation
        vol_term = 1.0 - nA * b_AA - nB * b_BB
        if vol_term <= 1e-9: # Density is at or beyond the physical limit
            return [np.inf, np.inf]

        # Calculate elements of the Jacobian matrix J_ij = ∂μ_i/∂n_j
        J11 = R * T / nA + R * T * b_AA**2 / vol_term**2 - 2 * a_AA
        J12 = R * T * b_AA * b_BB / vol_term**2 - 2 * a_AB
        J21 = R * T * b_BB * b_AA / vol_term**2 - 2 * a_AB
        J22 = R * T / nB + R * T * b_BB**2 / vol_term**2 - 2 * a_BB
        J = np.array([[J11, J12], [J21, J22]])

        # Gravitational force term vector
        F_g = -g * np.array([M_A, M_B])

        # Solve J * (dn/dz) = F_g for dn/dz
        try:
            dndt = np.linalg.solve(J, F_g)
            return dndt
        except np.linalg.LinAlgError:
            # Return large values if Jacobian is singular (e.g., at a critical point)
            return [np.inf, np.inf]

    # 3. Define the objective function for the root-finding algorithm
    def objective(n0):
        """
        Calculates the error between integrated total moles and actual total moles
        for a given set of initial densities n0 = [nA(0), nB(0)].
        """
        # Integrate the ODE system from z=0 to z=H
        sol = solve_ivp(get_derivatives, [0, H], n0, dense_output=True, method='RK45')

        # Return a large error if the integration failed
        if sol.status != 0:
            return [1e6, 1e6]

        # Evaluate the density profiles on a fine grid for integration
        z_eval = np.linspace(0, H, 200)
        n_profiles = sol.sol(z_eval)
        nA_profile, nB_profile = n_profiles[0, :], n_profiles[1, :]

        # Return a large error for non-physical negative densities
        if np.any(nA_profile < 0) or np.any(nB_profile < 0):
            return [1e6, 1e6]

        # Calculate total moles by integrating n(z)*A over the volume
        n_A_calc = trapz(nA_profile, z_eval) * A
        n_B_calc = trapz(nB_profile, z_eval) * A
        
        # Return the difference (error) to be minimized by the root finder
        return [n_A_calc - n_tot_A, n_B_calc - n_tot_B]

    # 4. Find the correct initial densities n(0) using the shooting method
    # Initial guess: average density, scaled up slightly as density is higher at the bottom
    V_total = A * H
    initial_guess = [
        1.2 * n_tot_A / V_total,
        1.2 * n_tot_B / V_total
    ]
    
    solution = root(objective, initial_guess, method='hybr', options={'xtol': 1e-8})
    if not solution.success:
        print("Warning: Root-finding for initial conditions may not have converged.")
        print("Solver message:", solution.message)
    
    n0_solution = solution.x

    # 5. Final calculation and output
    print("Density Profile Calculation Results")
    print("The final mass density profile is computed using the equation:")
    print(f"ρ(z) = n_A(z) * {M_A} + n_B(z) * {M_B}")
    print("-" * 75)
    print(f"{'Height z (m)':<15} | {'n_A(z) (mol/m^3)':<20} | {'n_B(z) (mol/m^3)':<20} | {'ρ(z) (kg/m^3)':<20}")
    print("-" * 75)
    
    # Generate profile points for final output table
    z_final = np.linspace(0, H, 11)
    final_sol = solve_ivp(get_derivatives, [0, H], n0_solution, t_eval=z_final)

    nA_z = final_sol.y[0]
    nB_z = final_sol.y[1]

    # Calculate the mass density profile rho(z)
    rho_z = nA_z * M_A + nB_z * M_B

    for i, z in enumerate(z_final):
        print(f"{z:<15.1f} | {nA_z[i]:<20.4f} | {nB_z[i]:<20.4f} | {rho_z[i]:<20.4f}")
        
    # For the final answer block, calculate the density at H/2
    rho_at_half_H = np.interp(H/2, z_final, rho_z)
    print(f"\n<<<The mass density at z = 5.0 m is {rho_at_half_H:.4f} kg/m^3.>>>")


if __name__ == '__main__':
    solve_density_profile()