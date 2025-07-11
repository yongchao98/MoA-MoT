import numpy as np
from scipy.optimize import brentq
from scipy.integrate import trapz

def solve_density_profile():
    """
    Calculates and prints the density profile of a non-ideal gas mixture
    in a vertical cylinder under gravity using a single-fluid hydrostatic model.
    """
    # --- 1. Define Constants and Parameters (in SI units) ---
    # Physical constants
    R = 8.31446  # Universal gas constant (J/(mol·K))
    N_Avo = 6.02214e23  # Avogadro's number (mol⁻¹)
    g = 9.81  # Gravitational acceleration (m/s²)

    # System parameters
    A = 0.1  # Cross-sectional area (m²)
    H = 10.0  # Height of the container (m)
    T = 500.0  # Temperature (K)

    # Gas A parameters
    N_A_particles = 2.0e23  # Number of particles of Gas A
    M_A = 28.0 / 1000  # Molar mass of Gas A (kg/mol)
    a_AA = 2.5  # vdW 'a' for Gas A (Pa·m⁶·mol⁻²)
    b_AA = 0.04  # vdW 'b' for Gas A (m³·mol⁻¹) - Assumed as given, though large

    # Gas B parameters
    N_B_particles = 1.5e23  # Number of particles of Gas B
    M_B = 44.0 / 1000  # Molar mass of Gas B (kg/mol)
    a_BB = 3.6  # vdW 'a' for Gas B (Pa·m⁶·mol⁻²)
    b_BB = 0.05  # vdW 'b' for Gas B (m³·mol⁻¹)

    # Interaction parameter
    a_AB = 3.0  # vdW 'a' for A-B interaction (Pa·m⁶·mol⁻²)

    # --- 2. Calculate Effective Mixture Properties ---
    N_total_particles = N_A_particles + N_B_particles
    N_total_moles = N_total_particles / N_Avo
    
    x_A = N_A_particles / N_total_particles
    x_B = N_B_particles / N_total_particles
    
    # Effective molar mass
    M_mix = x_A * M_A + x_B * M_B
    
    # Effective van der Waals parameters using mixing rules
    a_mix = x_A**2 * a_AA + 2 * x_A * x_B * a_AB + x_B**2 * a_BB
    b_mix = x_A * b_AA + x_B * b_BB
    
    # --- 3. Define the Governing Thermodynamic Function G(n) ---
    # This function comes from integrating the hydrostatic equation for a vdW gas
    def G_potential(n_m, a, b, temp):
        if n_m <= 0 or n_m >= 1/b:
            return np.nan # Return NaN for invalid densities
        term1 = np.log(n_m / (1 - n_m * b))
        term2 = 1 / (1 - n_m * b)
        term3 = (2 * a * n_m) / (temp * R)
        return temp * R * (term1 + term2) - 2 * a * n_m

    # --- 4. Numerical Solver for the Profile ---
    z_grid = np.linspace(0, H, 101)  # Discretize height into 100 segments

    def calculate_total_moles(n_m0):
        """
        Calculates the total moles in the cylinder for a given density at the bottom (n_m0).
        This is the function we want to find the root of.
        """
        n_m_profile = np.zeros_like(z_grid)
        
        # Check if n_m0 is valid
        if n_m0 <= 0 or n_m0 >= 1/b_mix:
            return np.inf

        # Value of potential at z=0
        G0 = G_potential(n_m0, a_mix, b_mix, T)
        n_m_profile[0] = n_m0
        
        # Solve for n_m(z) at each height z
        for i, z in enumerate(z_grid[1:], 1):
            target_G = G0 - M_mix * g * z
            
            # Function to find root of for n_m(z)
            def find_n_at_z(n):
                return G_potential(n, a_mix, b_mix, T) - target_G
            
            # Find n_m(z) using a bounded root solver
            try:
                # The density at the previous step is a good guess
                upper_bound = n_m_profile[i-1]
                n_m_profile[i] = brentq(find_n_at_z, 1e-9, upper_bound)
            except (ValueError, RuntimeError):
                # If root finding fails (e.g., no solution in the bracket)
                return np.inf # Return a large number to indicate failure

        # Integrate the profile to get total moles
        calculated_moles = trapz(n_m_profile, z_grid) * A
        return calculated_moles

    def objective_function(n_m0):
        return calculate_total_moles(n_m0) - N_total_moles

    # --- 5. Find the Correct Bottom Density n_m(0) ---
    # Find a search bracket for the root-finder
    n_avg = N_total_moles / (A * H)
    low_guess = n_avg * 0.1
    high_guess = n_avg * 10
    
    # Adjust bracket if objective function values have the same sign
    while np.sign(objective_function(low_guess)) == np.sign(objective_function(high_guess)):
        low_guess /= 2
        high_guess *= 2
        if high_guess >= 1/b_mix:
            high_guess = 0.999/b_mix
            break
            
    n_m0_solution = brentq(objective_function, low_guess, high_guess)

    # --- 6. Calculate Final Density Profile and Print Results ---
    final_n_m_profile = np.zeros_like(z_grid)
    G0_final = G_potential(n_m0_solution, a_mix, b_mix, T)
    final_n_m_profile[0] = n_m0_solution

    for i, z in enumerate(z_grid[1:], 1):
        target_G = G0_final - M_mix * g * z
        def find_n_at_z_final(n):
            return G_potential(n, a_mix, b_mix, T) - target_G
        final_n_m_profile[i] = brentq(find_n_at_z_final, 1e-9, final_n_m_profile[i-1])

    rho_profile = final_n_m_profile * M_mix

    # --- Final Output ---
    print("This script calculates the density profile ρ(z) of a non-ideal gas mixture under gravity.")
    print("The model uses hydrostatic equilibrium with the van der Waals equation for an effective single fluid.")
    print("\nThe final density is given by the equation:")
    print(f"  ρ(z) = n_m(z) * M_mix")
    print("\nWhere the molar density n_m(z) is solved numerically from the integrated hydrostatic equation:")
    print("  G(n_m(z)) - G(n_m(0)) = -M_mix * g * z")
    print("  with G(n) = RT[ln(n/(1-nb)) + 1/(1-nb)] - 2an")
    
    print("\n--- Calculated Effective Mixture Parameters ---")
    print(f"Total Moles: {N_total_moles:.4f} mol")
    print(f"M_mix (Effective Molar Mass): {M_mix:.5f} kg/mol")
    print(f"a_mix (Effective Attraction): {a_mix:.4f} Pa·m⁶·mol⁻²")
    print(f"b_mix (Effective Excluded Volume): {b_mix:.5f} m³·mol⁻¹")
    print(f"n_m(0) (Molar Density at z=0): {n_m0_solution:.4f} mol/m³")

    print("\n--- Mass Density Profile ρ(z) ---")
    print("Height (m) | Density (kg/m³)")
    print("-----------|----------------")
    for i in range(0, len(z_grid), 10):
        z_val = z_grid[i]
        rho_val = rho_profile[i]
        print(f"{z_val:^10.2f} | {rho_val:^15.5f}")

solve_density_profile()