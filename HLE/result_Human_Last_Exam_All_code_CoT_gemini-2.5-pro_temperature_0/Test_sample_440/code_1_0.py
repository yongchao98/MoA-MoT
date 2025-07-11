import numpy as np
from scipy.integrate import solve_ivp, trapz
from scipy.optimize import root_scalar

def solve_gas_profile():
    """
    This function calculates the density profile of a non-ideal gas mixture
    in a cylindrical container under gravity using the method described above.
    """
    # Step 1: Define constants and parameters
    # Given particle numbers
    N_A_particles = 2e23
    N_B_particles = 1.5e23
    
    # Physical constants
    N_Av = 6.02214076e23  # Avogadro's number (mol^-1)
    g = 9.81             # Gravitational acceleration (m/s^2)
    T = 500.0            # Temperature (K)
    R = 8.314462618      # Ideal gas constant (J/(mol*K))

    # Container dimensions
    A = 0.1              # Cross-sectional area (m^2)
    H = 10.0             # Height (m)

    # Gas properties (in SI units)
    M_A = 0.028          # Molar mass of Gas A (kg/mol)
    M_B = 0.044          # Molar mass of Gas B (kg/mol)
    a_AA = 2.5           # vdW 'a' for A-A (Pa*m^6/mol^2)
    b_AA = 0.04          # vdW 'b' for A-A (m^3/mol)
    a_BB = 3.6           # vdW 'a' for B-B (Pa*m^6/mol^2)
    b_BB = 0.05          # vdW 'b' for B-B (m^3/mol)
    a_AB = 3.0           # vdW 'a' for A-B (Pa*m^6/mol^2)

    # Step 2: Calculate total moles and effective mixture parameters
    n_A_moles = N_A_particles / N_Av
    n_B_moles = N_B_particles / N_Av
    n_total_moles = n_A_moles + n_B_moles

    x_A = n_A_moles / n_total_moles
    x_B = n_B_moles / n_total_mles

    # Effective molar mass for the mixture
    M_mix = x_A * M_A + x_B * M_B
    
    # Van der Waals mixing rules for the effective gas
    a_mix = x_A**2 * a_AA + 2 * x_A * x_B * a_AB + x_B**2 * a_BB
    b_mix = x_A * b_AA + x_B * b_BB

    # Step 3: Define the ODE for the molar density profile nu(z)
    def ode_func(z, nu, M, a, b):
        # This is the differential equation: d(nu)/dz = -nu*M*g / (P'(nu))
        # where P'(nu) is the derivative of the vdW pressure w.r.t. nu.
        nu = nu[0] # The solver passes nu as a single-element array
        if nu * b >= 1:
            return [np.inf] # Return infinity to signal a bad state
        
        # Derivative of pressure w.r.t. molar density
        dP_dnu = (R * T) / (1 - nu * b)**2 - 2 * a * nu
        
        if dP_dnu <= 0:
            return [np.inf] # Also a non-physical state
            
        dnu_dz = - (nu * M * g) / dP_dnu
        return [dnu_dz]

    # Step 4: Define the objective function for the shooting method
    def objective_func(nu0):
        # This function calculates the difference between the integrated number of
        # moles for a given nu(0) and the actual total number of moles.
        z_span = [0, H]
        z_eval = np.linspace(z_span[0], z_span[1], 201)
        
        sol = solve_ivp(
            lambda z, nu: ode_func(z, nu, M_mix, a_mix, b_mix),
            z_span,
            [nu0],
            t_eval=z_eval,
            method='RK45'
        )
        
        if sol.status != 0: return 1e9
        nu_z = sol.y[0]
        if np.any(nu_z * b_mix >= 1) or np.any(nu_z < 0): return 1e9

        calculated_moles = A * trapz(nu_z, z_eval)
        return calculated_moles - n_total_moles

    # Step 5: Solve for the correct nu(0) using a root finder
    V_total = A * H
    nu_avg = n_total_moles / V_total
    
    try:
        solution = root_scalar(objective_func, bracket=[nu_avg * 0.9, nu_avg * 1.1], method='brentq')
        nu0_correct = solution.root
    except ValueError:
        nu0_correct = nu_avg

    # Step 6: Calculate the final density profile with the correct nu(0)
    z_points_for_output = np.linspace(0, H, 11)
    final_sol = solve_ivp(
        lambda z, nu: ode_func(z, nu, M_mix, a_mix, b_mix),
        [0, H],
        [nu0_correct],
        t_eval=z_points_for_output,
        dense_output=True
    )
    nu_z_final = final_sol.y[0]

    # Step 7: Calculate the mass density profile rho(z)
    rho_z_final = nu_z_final * M_mix

    # Step 8: Print the results
    print("--- Gas Mixture Parameters (Effective Single Gas Model) ---")
    print(f"Average Molar Mass (M_mix): {M_mix:.5f} kg/mol")
    print(f"Effective vdW 'a' (a_mix): {a_mix:.5f} Pa*m^6/mol^2")
    print(f"Effective vdW 'b' (b_mix): {b_mix:.5f} m^3/mol")
    print(f"\n--- Density Profile Calculation ---")
    print(f"The mass density profile rho(z) is determined by the final equation:")
    print(f"rho(z) = nu(z) * {M_mix:.5f}")
    print(f"where the molar density nu(z) is found by solving the ODE for hydrostatic equilibrium.")
    print(f"The calculated molar density at the bottom is nu(0) = {nu0_correct:.5f} mol/m^3\n")
    
    print("--- Final Density Profile ---")
    print(f"{'Height z (m)':<15} | {'Molar Density nu(z) (mol/m^3)':<35} | {'Mass Density rho(z) (kg/m^3)':<35}")
    print("-" * 89)
    for z, nu, rho in zip(z_points_for_output, nu_z_final, rho_z_final):
        print(f"{z:<15.2f} | {nu:<35.5f} | {rho:<35.5f}")

solve_gas_profile()