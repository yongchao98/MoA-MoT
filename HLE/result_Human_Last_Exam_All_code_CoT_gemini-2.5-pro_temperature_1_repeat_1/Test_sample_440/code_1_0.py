import math

def solve_gas_density_profile():
    """
    Calculates and presents the density profile of a non-ideal gas mixture
    in a vertical cylinder under gravity.
    """
    # --- 1. Parameters and Constants ---
    # Given parameters
    A = 0.1  # Cross-sectional area, m^2
    H = 10.0  # Height of container, m
    T = 500.0  # Temperature, K
    g = 9.81  # Gravitational acceleration, m/s^2
    N_particles_A = 2e23
    N_particles_B = 1.5e23
    M_A_gmol = 28.0  # Molar mass of Gas A, g/mol
    M_B_gmol = 44.0  # Molar mass of Gas B, g/mol
    a_AA = 2.5  # van der Waals a, Pa*m^6/mol^2
    b_AA = 0.04  # van der Waals b, m^3/mol
    a_BB = 3.6
    b_BB = 0.05
    a_AB = 3.0

    # Physical constants
    N_AVOGADRO = 6.022e23  # mol^-1
    R = 8.314  # Gas constant, J/(mol*K)

    # Convert units to SI kg/mol
    M_A = M_A_gmol / 1000.0  # kg/mol
    M_B = M_B_gmol / 1000.0  # kg/mol

    # --- 2. Mixture Properties (Constant Composition Approximation) ---
    n_moles_A = N_particles_A / N_AVOGADRO
    n_moles_B = N_particles_B / N_AVOGADRO
    n_moles_total = n_moles_A + n_moles_B

    x_A = n_moles_A / n_moles_total
    x_B = n_moles_B / n_moles_total

    # Effective molar mass
    M_eff = x_A * M_A + x_B * M_B

    # Effective van der Waals parameters for the mixture
    a_mix = x_A**2 * a_AA + 2 * x_A * x_B * a_AB + x_B**2 * a_BB
    b_mix = x_A * b_AA + x_B * b_BB

    # --- 3. Numerical Method Implementation ---
    # Define the ODE: d(n_m)/dz = f(z, n_m)
    def get_dnm_dz(n_m):
        # This function calculates the derivative of molar density w.r.t. height
        # based on hydrostatic equilibrium and the van der Waals equation.
        denominator = (R * T / (1 - n_m * b_mix)**2) - (2 * a_mix * n_m)
        if denominator == 0:
            return float('inf')
        numerator = -n_m * M_eff * g
        return numerator / denominator

    # Define the objective function for the root finder
    def calculate_moles_error(n0_guess, dz=0.1):
        # Solves the ODE with a guess for n_m(0) and returns the error in total moles.
        z_vals = [0.0]
        n_m_vals = [n0_guess]
        
        # Euler method to integrate the ODE
        current_z = 0.0
        while current_z < H:
            n_m = n_m_vals[-1]
            if n_m <= 0 or (1 - n_m * b_mix) <= 0: # Physicality check
                return float('inf') # Return a large error if density becomes non-physical
            dnm_dz = get_dnm_dz(n_m)
            n_m_new = n_m + dnm_dz * dz
            current_z += dz
            n_m_vals.append(n_m_new)
            z_vals.append(current_z)

        # Trapezoidal rule to integrate the profile for total moles
        integrated_moles = 0
        for i in range(len(z_vals) - 1):
            integrated_moles += (n_m_vals[i] + n_m_vals[i+1]) / 2.0 * (z_vals[i+1] - z_vals[i])
        
        calculated_total_moles = A * integrated_moles
        return calculated_total_moles - n_moles_total

    # Bisection method to find the root of the objective function
    # This finds the correct initial density n_m(0)
    low_guess = n_moles_total / (A * H) * 0.1 # Start with a low guess
    high_guess = n_moles_total / (A * H) * 10.0 # and a high guess
    
    # Ensure the initial bracket is valid
    while calculate_moles_error(high_guess) < 0:
        high_guess *= 2
    while calculate_moles_error(low_guess) > 0:
        low_guess /= 2

    for _ in range(100): # 100 iterations for high precision
        mid_guess = (low_guess + high_guess) / 2
        error = calculate_moles_error(mid_guess)
        if abs(error) < 1e-9:
            break
        if error < 0:
            low_guess = mid_guess
        else:
            high_guess = mid_guess
    
    n_m_0 = (low_guess + high_guess) / 2.0

    # --- 4. Final Profile Calculation and Output ---
    # Solve the ODE one last time with the correct initial condition
    z_profile = []
    n_m_profile = []
    dz_final = 0.05
    z, n_m = 0.0, n_m_0
    while z <= H:
        z_profile.append(z)
        n_m_profile.append(n_m)
        dnm_dz = get_dnm_dz(n_m)
        n_m += dnm_dz * dz_final
        z += dz_final

    # Calculate mass density profile
    rho_profile = [n * M_eff for n in n_m_profile]

    # Find densities at specific heights
    def get_rho_at_z(z_target):
        # Simple interpolation/nearest neighbor
        idx = min(range(len(z_profile)), key=lambda i: abs(z_profile[i] - z_target))
        return rho_profile[idx]

    rho_0 = get_rho_at_z(0.0)
    rho_mid = get_rho_at_z(H / 2.0)
    rho_H = get_rho_at_z(H)

    # Print the results in the required format
    print("The mass density profile rho(z) is determined by numerically solving the hydrostatic equilibrium equation dP/dz = -rho(z)*g, using the van der Waals equation for a gas mixture under a constant composition approximation.")
    print("\nThe governing differential equation for the molar density n_m(z) is:")
    print("d(n_m)/dz = - (n_m * M_eff * g) / (R * T / (1 - n_m * b_mix)^2 - 2 * a_mix * n_m)")
    
    print("\nThe parameters for this equation are:")
    print(f"Effective Molar Mass (M_eff): {M_eff:.5f} kg/mol")
    print(f"Gravitational Acceleration (g): {g:.2f} m/s^2")
    print(f"Gas Constant (R): {R:.3f} J/(mol*K)")
    print(f"Temperature (T): {T:.1f} K")
    print(f"Mixture Attraction Parameter (a_mix): {a_mix:.5f} Pa*m^6/mol^2")
    print(f"Mixture Excluded Volume (b_mix): {b_mix:.5f} m^3/mol")

    print("\nThe initial condition n_m(z=0) is found by requiring the total number of moles in the container to be constant.")
    print(f"Known Total Moles (N_total): {n_moles_total:.5f} mol")
    
    print("\nA numerical shooting method yields the required initial molar density:")
    print(f"n_m(z=0) = {n_m_0:.5f} mol/m^3")
    
    print("\nFinally, the mass density profile, rho(z) = n_m(z) * M_eff, is calculated. At selected heights, the density is:")
    print(f"rho(z = {0.0:.1f} m) = {rho_0:.5f} kg/m^3")
    print(f"rho(z = {H/2.0:.1f} m) = {rho_mid:.5f} kg/m^3")
    print(f"rho(z = {H:.1f} m) = {rho_H:.5f} kg/m^3")


solve_gas_density_profile()