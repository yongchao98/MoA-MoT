import numpy as np
from scipy.integrate import solve_ivp, trapz
from scipy.optimize import root

def solve_density_profile():
    """
    Solves for the density profile of a two-component van der Waals gas
    mixture in a gravitational field.
    """
    # --- 1. Parameters and Constants ---
    # Given parameters
    params = {
        'A': 0.1,  # Cross-sectional area (m^2)
        'H': 10.0,  # Height of container (m)
        'T': 500.0,  # Temperature (K)
        'g': 9.81,  # Gravitational acceleration (m/s^2)
        'N_A': 2e23,  # Number of particles of Gas A
        'N_B': 1.5e23,  # Number of particles of Gas B
        'M_A': 28e-3,  # Molar mass of Gas A (kg/mol)
        'M_B': 44e-3,  # Molar mass of Gas B (kg/mol)
        'a_AA': 2.5,  # vdW a for A-A (Pa m^6 mol^-2)
        'b_AA': 0.04 / 1000, # vdW b for A-A (m^3 mol^-1), converted from L/mol
        'a_BB': 3.6,  # vdW a for B-B (Pa m^6 mol^-2)
        'b_BB': 0.05 / 1000, # vdW b for B-B (m^3 mol^-1), converted from L/mol
        'a_AB': 3.0,  # vdW a for A-B (Pa m^6 mol^-2)
        # Constants
        'R': 8.314,  # Ideal gas constant (J mol^-1 K^-1)
        'N_avogadro': 6.022e23  # Avogadro's number (mol^-1)
    }

    # Total moles of each gas
    n_A_total = params['N_A'] / params['N_avogadro']
    n_B_total = params['N_B'] / params['N_avogadro']
    params['n_A_total'] = n_A_total
    params['n_B_total'] = n_B_total

    # --- 2. Define the ODE System ---
    # The system is derived from the condition of hydrostatic and diffusive equilibrium:
    # d(mu_i)/dz = -M_i * g
    # This leads to a matrix equation: J * [dn_A/dz, dn_B/dz]^T = [-M_A*g, -M_B*g]^T
    # where J is the Jacobian matrix J_ij = d(mu_i)/d(n_j)

    def ode_system(z, y, p):
        n_A, n_B = y
        
        # Unpack parameters for clarity
        R, T = p['R'], p['T']
        M_A, M_B = p['M_A'], p['M_B']
        g = p['g']
        a_AA, a_BB, a_AB = p['a_AA'], p['a_BB'], p['a_AB']
        b_AA, b_BB = p['b_AA'], p['b_BB']

        # Check for unphysical densities
        if n_A <= 0 or n_B <= 0:
            return [0, 0]

        n_t = n_A + n_B
        b_lin = n_A * b_AA + n_B * b_BB
        D = 1.0 - b_lin

        # Check if the density is too high (volume of molecules exceeds container volume)
        if D <= 1e-6: # Avoid division by zero
            return [0, 0]

        # Jacobian elements J_ij = d(mu_i)/d(n_j)
        # Based on mu_i = mu_i^0 + RT*ln(n_i) - RT*ln(1-b_lin) + (n_t*RT*b_i)/(1-b_lin) - 2*sum(n_j*a_ij)
        J_AA = R*T/n_A + 2*R*T*b_AA/D + R*T*n_t*b_AA**2/D**2 - 2*a_AA
        J_BB = R*T/n_B + 2*R*T*b_BB/D + R*T*n_t*b_BB**2/D**2 - 2*a_BB
        J_AB = R*T*(b_AA+b_BB)/D + R*T*n_t*b_AA*b_BB/D**2 - 2*a_AB
        J_BA = J_AB # Jacobian must be symmetric

        # Solve the 2x2 linear system for [dnA/dz, dnB/dz]
        det_J = J_AA * J_BB - J_AB * J_BA
        if abs(det_J) < 1e-9: # Avoid singularity (could indicate phase transition)
             return [0,0]

        # Gravitational force terms
        F_A = -M_A * g
        F_B = -M_B * g

        # Cramer's rule
        dnA_dz = (F_A * J_BB - F_B * J_AB) / det_J
        dnB_dz = (F_B * J_AA - F_A * J_AB) / det_J

        return [dnA_dz, dnB_dz]

    # --- 3. Objective Function for Root Finding (Shooting Method) ---
    def objective_func(y0, p):
        n_A0, n_B0 = y0
        if n_A0 < 0 or n_B0 < 0:
            return [1e6, 1e6] # Large error for invalid initial guess

        sol = solve_ivp(ode_system, [0, p['H']], y0, args=(p,), dense_output=True, method='RK45')
        
        # Check if integration was successful
        if sol.status != 0:
            return [1e6, 1e6]

        z_eval = np.linspace(0, p['H'], 200)
        n_profiles = sol.sol(z_eval)
        
        n_A_calc = trapz(n_profiles[0], z_eval) * p['A']
        n_B_calc = trapz(n_profiles[1], z_eval) * p['A']
        
        err_A = n_A_calc - p['n_A_total']
        err_B = n_B_calc - p['n_B_total']
        
        return [err_A, err_B]

    # --- 4. Solve the Boundary Value Problem ---
    # Initial guess for molar densities at z=0
    V_total = params['A'] * params['H']
    n_A_avg = n_A_total / V_total
    n_B_avg = n_B_total / V_total
    # Density at the bottom should be higher than average
    y0_guess = [1.5 * n_A_avg, 1.5 * n_B_avg]

    # Find the correct initial densities n_A(0), n_B(0)
    root_sol = root(objective_func, y0_guess, args=(params,), method='hybr')
    
    if not root_sol.success:
        print("Root finding failed to converge. The result might be inaccurate.")
        print(f"Message: {root_sol.message}")
        return

    y0_final = root_sol.x

    # --- 5. Final Calculation and Output ---
    # Integrate one last time with the correct initial conditions
    final_sol = solve_ivp(ode_system, [0, params['H']], y0_final, args=(params,), dense_output=True, num_points=101)

    z_points = final_sol.t
    n_A_z, n_B_z = final_sol.y

    # Calculate the total mass density profile rho(z)
    rho_z = n_A_z * params['M_A'] + n_B_z * params['M_B']

    # Print the governing equations and results
    print("The density profile is found by solving the hydrostatic and diffusive equilibrium equations:")
    print("d(mu_A)/dz = -M_A * g")
    print("d(mu_B)/dz = -M_B * g")
    print("\nWhich can be written in matrix form as:")
    print("[[J_AA, J_AB], [J_BA, J_BB]] * [dn_A/dz, dn_B/dz]^T = [-M_A*g, -M_B*g]^T")
    print("\nFor this system, the numerical values on the right-hand side are:")
    print(f"-M_A * g = -{params['M_A']:.3f} kg/mol * {params['g']:.2f} m/s^2 = {(-params['M_A'] * params['g']):.4f} N/mol")
    print(f"-M_B * g = -{params['M_B']:.3f} kg/mol * {params['g']:.2f} m/s^2 = {(-params['M_B'] * params['g']):.4f} N/mol")
    print("\nThe Jacobian J depends on local densities n_A(z) and n_B(z) via the van der Waals model.")

    print("\n--- Numerical Solution for the Mass Density Profile rho(z) ---")
    
    # Find indices for z=0, H/2, H
    idx_0 = np.argmin(np.abs(z_points - 0.0))
    idx_mid = np.argmin(np.abs(z_points - params['H']/2.0))
    idx_H = np.argmin(np.abs(z_points - params['H']))
    
    rho_0 = rho_z[idx_0]
    rho_mid = rho_z[idx_mid]
    rho_H = rho_z[idx_H]

    print(f"The equation for the mass density profile rho(z) is numerically solved.")
    print("Here are its values at key heights:")
    print(f"rho(z = {z_points[idx_0]:.1f} m) = {rho_0:.4f} kg/m^3")
    print(f"rho(z = {z_points[idx_mid]:.1f} m) = {rho_mid:.4f} kg/m^3")
    print(f"rho(z = {z_points[idx_H]:.1f} m) = {rho_H:.4f} kg/m^3")
    
    # Final answer format
    print(f"\n<<<rho(z = 5.0 m) = {rho_mid:.4f} kg/m^3>>>")

if __name__ == '__main__':
    solve_density_profile()