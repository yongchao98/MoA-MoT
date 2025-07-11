import numpy as np
from scipy.integrate import solve_ivp, trapz
from scipy.optimize import root

def solve_density_profile():
    """
    Calculates and prints the density profile of a two-component non-ideal gas
    mixture in a gravitational field using the van der Waals model.
    """
    # --- 1. Parameters and Constants ---
    # Given parameters in SI units
    A = 0.1  # Cross-sectional area (m^2)
    H = 10.0  # Height of container (m)
    T = 500.0  # Temperature (K)
    g = 9.81  # Gravitational acceleration (m/s^2)
    
    # Gas A properties
    N_A = 2.0e23  # Number of particles
    M_A = 0.028  # Molar mass (kg/mol)
    a_AA = 2.5  # van der Waals 'a' (Pa * m^6 * mol^-2)
    b_AA = 0.04  # van der Waals 'b' (m^3 * mol^-1)

    # Gas B properties
    N_B = 1.5e23  # Number of particles
    M_B = 0.044  # Molar mass (kg/mol)
    a_BB = 3.6  # van der Waals 'a' (Pa * m^6 * mol^-2)
    b_BB = 0.05  # van der Waals 'b' (m^3 * mol^-1)

    # Interaction parameter
    a_AB = 3.0  # (Pa * m^6 * mol^-2)

    # Physical constants
    R = 8.31446  # Ideal gas constant (J * mol^-1 * K^-1)
    N_AVOGADRO = 6.02214e23  # Avogadro's number (mol^-1)

    # Target total moles for each gas
    N_A_tot_target = N_A / N_AVOGADRO
    N_B_tot_target = N_B / N_AVOGADRO

    # --- 2. Define the System of ODEs ---
    # Based on d(mu_i)/dz = -M_i * g, where mu is the chemical potential.
    # This function calculates [dc_A/dz, dc_B/dz]
    def odes(z, c):
        c_A, c_B = c

        # Concentrations must be positive. Return large gradients to steer solver away.
        if c_A <= 1e-6 or c_B <= 1e-6:
            return [1e6, 1e6]

        # Denominator in vdW terms, represents free volume factor
        denom_b = 1.0 - c_A * b_AA - c_B * b_BB
        
        # Unphysical state (negative free volume).
        if denom_b <= 1e-6:
            return [1e6, 1e6]
            
        V_free_inv = 1.0 / denom_b

        # Calculate the Jacobian matrix J_ij = d(mu_i)/d(c_j)
        # These are the partial derivatives of the chemical potentials
        J_AA = (R * T / c_A) + (R * T * b_AA * V_free_inv) + (R * T * b_AA**2 * c_A * V_free_inv**2) - 2 * a_AA
        J_AB = (R * T * b_BB * V_free_inv) + (R * T * b_AA * b_BB * c_A * V_free_inv**2) - 2 * a_AB
        J_BA = (R * T * b_AA * V_free_inv) + (R * T * b_BB * b_AA * c_B * V_free_inv**2) - 2 * a_AB
        J_BB = (R * T / c_B) + (R * T * b_BB * V_free_inv) + (R * T * b_BB**2 * c_B * V_free_inv**2) - 2 * a_BB

        # Calculate the determinant of the Jacobian
        det_J = J_AA * J_BB - J_AB * J_BA
        if abs(det_J) < 1e-9: # Avoid singularity
            return [0, 0]

        # Solve the linear system J * (dc/dz) = -M*g for dc/dz
        dc_A_dz = (-g * M_A * J_BB + g * M_B * J_AB) / det_J
        dc_B_dz = (g * M_A * J_BA - g * M_B * J_AA) / det_J

        return [dc_A_dz, dc_B_dz]

    # --- 3. Objective Function for Root Finding (Shooting Method) ---
    # This function returns the error in total moles for a given set of initial concentrations.
    def objective(c0):
        c_A0, c_B0 = c0
        if c_A0 <= 0 or c_B0 <= 0:
            return [1e6, 1e6] # Return large error for non-physical guess

        # Solve the initial value problem
        sol = solve_ivp(odes, [0, H], [c_A0, c_B0], dense_output=True, method='Radau')

        if not sol.success:
            return [1e6, 1e6] # Return large error if integration fails

        # Integrate the resulting profiles to find the total number of moles
        z_pts = np.linspace(0, H, 200)
        c_profiles = sol.sol(z_pts)
        
        # Check for unphysical results during integration
        if np.any(c_profiles < 0):
            return [1e6, 1e6]

        N_A_calc = trapz(c_profiles[0], z_pts) * A
        N_B_calc = trapz(c_profiles[1], z_pts) * A

        # Return the difference between calculated and target moles
        return [N_A_calc - N_A_tot_target, N_B_calc - N_B_tot_target]

    # --- 4. Find the Correct Initial Concentrations ---
    # Initial guess based on average molar density
    V_tot = A * H
    c_A_avg = N_A_tot_target / V_tot
    c_B_avg = N_B_tot_target / V_tot
    initial_guess = [c_A_avg, c_B_avg]

    # Use a root finder to find c0 that makes the objective function zero
    solution = root(objective, initial_guess, method='hybr', options={'xtol': 1e-8})
    if not solution.success:
        print("Warning: Root finding for initial conditions may not have converged.")
        print(solution.message)
    
    c_A0_final, c_B0_final = solution.x

    # --- 5. Calculate and Print the Final Density Profile ---
    # Solve the ODEs one last time with the correct initial conditions
    z_eval = np.linspace(0, H, 11) # Evaluate at 11 points including 0 and H
    final_sol = solve_ivp(odes, [0, H], [c_A0_final, c_B0_final], dense_output=True, t_eval=z_eval)

    c_A_final = final_sol.y[0]
    c_B_final = final_sol.y[1]

    # Calculate the mass density profile rho(z) = c_A(z)*M_A + c_B(z)*M_B
    rho_z = c_A_final * M_A + c_B_final * M_B

    print("Final Density Profile rho(z):")
    print("-" * 35)
    print(f"{'Height z (m)':<15} | {'Density rho (kg/m^3)':<20}")
    print("-" * 35)
    for z, rho in zip(final_sol.t, rho_z):
        print(f"{z:<15.1f} | {rho:<20.4f}")
    print("-" * 35)

    # Print the breakdown of the calculation for z=0 as requested
    print("\nExample: Calculation for z = 0.0 m")
    print(f"The density equation is: rho(z) = c_A(z) * M_A + c_B(z) * M_B")
    print(f"At z = 0.0 m:")
    print(f"rho(0) = {c_A_final[0]:.4f} [mol/m^3] * {M_A:.4f} [kg/mol] + {c_B_final[0]:.4f} [mol/m^3] * {M_B:.4f} [kg/mol]")
    rho_A_0 = c_A_final[0] * M_A
    rho_B_0 = c_B_final[0] * M_B
    print(f"rho(0) = {rho_A_0:.4f} [kg/m^3] + {rho_B_0:.4f} [kg/m^3]")
    print(f"rho(0) = {rho_z[0]:.4f} [kg/m^3]")


if __name__ == '__main__':
    solve_density_profile()