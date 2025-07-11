import numpy as np
from scipy.optimize import root, fsolve
from scipy.integrate import trapz

def solve_density_profile():
    """
    Calculates the density profile of a non-ideal gas mixture in a gravitational field.
    """
    # --- Parameters in SI units ---
    A = 0.1  # Cross-sectional area (m^2)
    H = 10.0  # Height of container (m)
    T = 500.0  # Temperature (K)
    g = 9.81  # Gravitational acceleration (m/s^2)
    
    # Gas A
    M_A = 28e-3  # Molar mass (kg/mol)
    N_A = 2e23   # Number of particles
    a_AA = 2.5   # Pa * m^6 * mol^-2
    b_AA = 0.04e-3  # m^3 * mol^-1 (Note: converted problem spec from L/mol to m^3/mol assuming typo)
                    # Let's stick to the problem statement:
    b_AA = 0.04 # m^3 * mol^-1, assuming units are already SI-consistent for V_m
    
    # Gas B
    M_B = 44e-3  # Molar mass (kg/mol)
    N_B = 1.5e23 # Number of particles
    a_BB = 3.6   # Pa * m^6 * mol^-2
    b_BB = 0.05  # m^3 * mol^-1

    # Interaction
    a_AB = 3.0  # Pa * m^6 * mol^-2

    # Physical constants
    R = 8.31446  # Ideal gas constant (J/mol*K)
    N_AVOGADRO = 6.02214e23  # Avogadro's number (mol^-1)
    
    # Total moles of each gas
    n_A_total = N_A / N_AVOGADRO
    n_B_total = N_B / N_AVOGADRO

    # --- Numerical simulation setup ---
    N_z = 201  # Number of points for discretization
    z_array = np.linspace(0, H, N_z)
    dz = H / (N_z - 1)

    def vdw_pressure(c_A, c_B):
        """Calculates VdW pressure for a mixture."""
        if c_A + c_B <= 0: return 0.0
        c_total = c_A + c_B
        x_A = c_A / c_total
        x_B = c_B / c_total

        # Mixing rules
        a_mix = x_A**2 * a_AA + x_B**2 * a_BB + 2 * x_A * x_B * a_AB
        b_mix = x_A * b_AA + x_B * b_BB
        
        # Van der Waals equation for molar concentrations c_A, c_B
        pressure_term = R * T * c_total / (1 - b_mix * c_total)
        attraction_term = a_mix * c_total**2
        
        if pressure_term < attraction_term:
             return 0.0 # avoid negative pressure
        return pressure_term - attraction_term

    def get_concentrations_from_pressure(P, z, K):
        """Finds c_B (and c_A) for a given pressure P, height z, and ratio K."""
        # Relative distribution due to gravity
        ratio_exp_term = np.exp(-(M_A - M_B) * g * z / (R * T))
        
        # Define the objective function for the root finder
        def pressure_error(c_B_guess):
            c_B_val = c_B_guess[0]
            if c_B_val < 0: return 1e9 # physically impossible
            c_A_val = c_B_val * K * ratio_exp_term
            return P - vdw_pressure(c_A_val, c_B_val)

        # Initial guess for c_B based on ideal gas law
        c_B_initial_guess = P / (R * T * (K * ratio_exp_term + 1))
        
        # Solve for c_B
        c_B_sol = fsolve(pressure_error, [c_B_initial_guess])
        c_B = c_B_sol[0]
        c_A = c_B * K * ratio_exp_term
        return c_A, c_B

    def objective_function(guess):
        """
        Calculates the error in total moles for a given P(0) and K.
        This function will be driven to zero by the root solver.
        """
        P0, K = guess
        if P0 < 0 or K < 0: return [1e9, 1e9]
        
        P = np.zeros(N_z)
        c_A = np.zeros(N_z)
        c_B = np.zeros(N_z)
        rho = np.zeros(N_z)

        P[0] = P0
        try:
            c_A[0], c_B[0] = get_concentrations_from_pressure(P[0], z_array[0], K)
            rho[0] = c_A[0] * M_A + c_B[0] * M_B
        except (ValueError, TypeError):
             return [1e9, 1e9]

        # Integrate up the cylinder
        for i in range(N_z - 1):
            # Update pressure using hydrostatic equilibrium
            P[i+1] = P[i] - rho[i] * g * dz
            if P[i+1] < 0: P[i+1] = 0

            # Find concentrations at the new pressure
            try:
                c_A[i+1], c_B[i+1] = get_concentrations_from_pressure(P[i+1], z_array[i+1], K)
                rho[i+1] = c_A[i+1] * M_A + c_B[i+1] * M_B
            except (ValueError, TypeError):
                # if solver fails, return large error
                return [1e9, 1e9]
        
        # Calculate total moles from the integrated profiles
        n_A_calc = trapz(c_A, z_array) * A
        n_B_calc = trapz(c_B, z_array) * A
        
        # Return the error vector
        return [n_A_calc - n_A_total, n_B_calc - n_B_total]

    # --- Solve the system ---
    # Initial guesses for P(0) and K
    avg_moles = (n_A_total + n_B_total)
    avg_vol = A*H
    P_guess = avg_moles * R * T / avg_vol
    K_guess = n_A_total / n_B_total
    initial_guess = [P_guess, K_guess]
    
    solution = root(objective_function, initial_guess, method='hybr')
    if not solution.success:
        print("Warning: Root finding did not converge. The results might be inaccurate.")
        print("Solver message:", solution.message)
    
    P0_sol, K_sol = solution.x

    # --- Generate the final profile with the correct parameters ---
    P_final = np.zeros(N_z)
    c_A_final = np.zeros(N_z)
    c_B_final = np.zeros(N_z)
    rho_final = np.zeros(N_z)

    P_final[0] = P0_sol
    c_A_final[0], c_B_final[0] = get_concentrations_from_pressure(P_final[0], z_array[0], K_sol)
    rho_final[0] = c_A_final[0] * M_A + c_B_final[0] * M_B

    for i in range(N_z - 1):
        P_final[i+1] = P_final[i] - rho_final[i] * g * dz
        c_A_final[i+1], c_B_final[i+1] = get_concentrations_from_pressure(P_final[i+1], z_array[i+1], K_sol)
        rho_final[i+1] = c_A_final[i+1] * M_A + c_B_final[i+1] * M_B

    # --- Print the final results ---
    print("Final Density Profile rho(z):")
    print("="*30)
    print(f"{'Height z (m)':<15} | {'Density rho (kg/m^3)':<20}")
    print("-"*30)
    
    # Print the "equation" by showing its value at several points
    for i in range(0, N_z, N_z // 10):
      z_val = z_array[i]
      rho_val = rho_final[i]
      print(f"At z = {z_val:5.2f},      rho = {rho_val:.4f}")
    
    # Make sure to print the last point
    if (N_z-1) % (N_z // 10) != 0:
        z_val = z_array[-1]
        rho_val = rho_final[-1]
        print(f"At z = {z_val:5.2f},      rho = {rho_val:.4f}")
    
    # To conform to the output format, provide a representative value
    # The density at the bottom of the cylinder, rho(0), is a good choice.
    print(f"\n<<<Density at z=0 is {rho_final[0]:.4f} kg/m^3>>>")

solve_density_profile()