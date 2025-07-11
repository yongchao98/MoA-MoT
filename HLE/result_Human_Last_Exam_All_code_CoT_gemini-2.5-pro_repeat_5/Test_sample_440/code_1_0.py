import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root

def solve_gas_profile():
    """
    Calculates and prints the density profile of a non-ideal gas mixture
    in a gravitational field.
    """
    # 1. Constants and Parameters (in SI units)
    # -----------------------------------------------------------------
    # System parameters
    A = 0.1  # Cross-sectional area (m^2)
    H = 10.0  # Height (m)
    T = 500.0  # Temperature (K)
    g = 9.81  # Gravitational acceleration (m/s^2)

    # Gas properties
    N_A_particles = 2.0e23
    N_B_particles = 1.5e23
    M_A_gmol = 28.0  # Molar mass of Gas A (g/mol)
    M_B_gmol = 44.0  # Molar mass of Gas B (g/mol)

    # Van der Waals parameters
    a_AA = 2.5  # Pa * m^6 * mol^-2
    a_BB = 3.6  # Pa * m^6 * mol^-2
    a_AB = 3.0  # Pa * m^6 * mol^-2
    b_A = 0.04  # m^3 * mol^-1 (Note: This is an unusually large value)
    b_B = 0.05  # m^3 * mol^-1 (Note: This is an unusually large value)

    # Universal constants
    R = 8.31446  # Universal gas constant (J * mol^-1 * K^-1)
    N_av = 6.02214e23  # Avogadro's number (mol^-1)

    # --- Convert to consistent SI units and create vectors/matrices ---
    M_A_kg = M_A_gmol / 1000.0  # kg/mol
    M_B_kg = M_B_gmol / 1000.0  # kg/mol
    M_vec = np.array([M_A_kg, M_B_kg])

    n_A_total = N_A_particles / N_av  # Total moles of A
    n_B_total = N_B_particles / N_av  # Total moles of B
    n_total_target = np.array([n_A_total, n_B_total])

    a_mat = np.array([[a_AA, a_AB], [a_AB, a_BB]])
    b_vec = np.array([b_A, b_B])

    # 2. Define the ODE System from dμ/dz = 0
    # -----------------------------------------------------------------
    def ode_system(z, n, R_const, T_const, g_const, M, a, b):
        """
        Defines the system of ODEs: d(n_i)/dz = f(n_A, n_B).
        n = [n_A, n_B] in mol/m^3
        """
        nA, nB = n
        if nA <= 0 or nB <= 0: return [np.inf, np.inf]

        D = 1 - (nA * b[0] + nB * b[1])
        if D <= 1e-6: return [np.inf, np.inf] # Unphysical density

        # Jacobian matrix J_ij = ∂μ_i/∂n_j
        J = np.zeros((2, 2))
        RT_D = R_const * T_const / D
        RT_D2 = R_const * T_const / (D**2)

        J[0, 0] = R_const * T_const / nA + RT_D * b[0] + RT_D2 * b[0]**2 - 2 * a[0, 0]
        J[0, 1] = RT_D * b[1] + RT_D2 * b[0] * b[1] - 2 * a[0, 1]
        J[1, 0] = RT_D * b[0] + RT_D2 * b[1] * b[0] - 2 * a[1, 0]
        J[1, 1] = R_const * T_const / nB + RT_D * b[1] + RT_D2 * b[1]**2 - 2 * a[1, 1]

        if np.linalg.det(J) == 0: return [np.inf, np.inf]
        
        J_inv = np.linalg.inv(J)
        dndt = -g_const * J_inv @ M
        return dndt

    # 3. Define the Objective Function for the Shooting Method
    # -----------------------------------------------------------------
    def objective_func(n0, R_const, T_const, g_const, M, a, b, A_cyl, H_cyl, n_target):
        """
        Calculates the error between integrated moles and target moles.
        """
        ode_func = lambda z, n: ode_system(z, n, R_const, T_const, g_const, M, a, b)
        sol = solve_ivp(ode_func, [0, H_cyl], n0, dense_output=True, method='RK45')

        if sol.status != 0 or np.any(sol.y < 0):
            return [1e10, 1e10] # Return large error if integration fails

        # Integrate profiles to get total calculated moles
        n_A_calc = A_cyl * np.trapz(sol.y[0], sol.t)
        n_B_calc = A_cyl * np.trapz(sol.y[1], sol.t)
        
        return [n_A_calc - n_target[0], n_B_calc - n_target[1]]

    # 4. Solve for the initial conditions n(0)
    # -----------------------------------------------------------------
    # Initial guess based on average density
    V_cyl = A * H
    n_avg_A = n_A_total / V_cyl
    n_avg_B = n_B_total / V_cyl
    # A slightly higher density at the bottom is a good starting guess
    initial_guess = [1.5 * n_avg_A, 1.5 * n_avg_B]
    
    # Wrap objective function to pass constant parameters
    obj_wrapper = lambda n0: objective_func(n0, R, T, g, M_vec, a_mat, b_vec, A, H, n_total_target)
    
    solution = root(obj_wrapper, initial_guess, method='hybr', tol=1e-8)
    
    if not solution.success:
        print("Root finding for initial conditions failed:", solution.message)
        return

    n0_solution = solution.x

    # 5. Calculate and Display the Final Density Profile
    # -----------------------------------------------------------------
    # Integrate one last time with the correct n(0) for high-resolution output
    final_ode_func = lambda z, n: ode_system(z, n, R, T, g, M_vec, a_mat, b_vec)
    z_eval_points = np.linspace(0, H, 11)
    final_sol = solve_ivp(final_ode_func, [0, H], n0_solution, dense_output=True, t_eval=z_eval_points)

    z_profile = final_sol.t
    n_A_profile = final_sol.y[0]
    n_B_profile = final_sol.y[1]
    
    # Calculate mass density profile in kg/m^3
    rho_profile = n_A_profile * M_A_kg + n_B_profile * M_B_kg

    print("--- Density Profile rho(z) ---")
    for z, rho in zip(z_profile, rho_profile):
        print(f"At z = {z:4.1f} m, Mass Density = {rho:.4f} kg/m^3")
    
    print("\n--- Detailed Calculation for the Final Equation at z=0 ---")
    rho_0 = rho_profile[0]
    nA_0 = n_A_profile[0]
    nB_0 = n_B_profile[0]
    
    print("The mass density is the sum of the partial mass densities:")
    print("rho(z) = n_A(z) * M_A + n_B(z) * M_B")
    print("\nAt z = 0.0 m:")
    print(f"rho(0) = n_A(0) * M_A + n_B(0) * M_B")
    print(f"rho(0) = ({nA_0:.5f} mol/m^3 * {M_A_kg:.3f} kg/mol) + ({nB_0:.5f} mol/m^3 * {M_B_kg:.3f} kg/mol)")
    print(f"rho(0) = {nA_0 * M_A_kg:.5f} kg/m^3 (from Gas A) + {nB_0 * M_B_kg:.5f} kg/m^3 (from Gas B)")
    print(f"rho(0) = {rho_0:.5f} kg/m^3")
    
    # The question is ambiguous about the final format, so we output the result at z=H/2
    rho_mid = np.interp(H/2, z_profile, rho_profile)
    print(f"\n<<<Density at H/2 (z=5.0m) is {rho_mid:.5f} kg/m^3>>>")

# Run the solver
solve_gas_profile()