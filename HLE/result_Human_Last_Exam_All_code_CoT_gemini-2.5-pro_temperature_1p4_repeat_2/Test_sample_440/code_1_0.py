import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.integrate import simps

def solve_density_profile():
    """
    Calculates and prints the density profile of a non-ideal gas mixture
    in a cylinder under gravity using the van der Waals equation and
    hydrostatic equilibrium.
    """
    # 1. Define Parameters and Constants
    # System Parameters
    A = 0.1         # Cross-sectional area (m^2)
    H = 10.0        # Height of container (m)
    T = 500.0       # Temperature (K)
    g = 9.81        # Gravitational acceleration (m/s^2)
    N_A = 2e23      # Number of particles of Gas A
    N_B = 1.5e23    # Number of particles of Gas B
    M_A_gmol = 28.0 # Molar mass of Gas A (g/mol)
    M_B_gmol = 44.0 # Molar mass of Gas B (g/mol)
    
    # van der Waals parameters
    a_AA = 2.5      # Pa * m^6 * mol^-2
    b_AA = 0.04     # m^3 * mol^-1
    a_BB = 3.6      # Pa * m^6 * mol^-2
    b_BB = 0.05     # m^3 * mol^-1
    a_AB = 3.0      # Pa * m^6 * mol^-2

    # Physical Constants
    R = 8.31446     # Universal gas constant (J/(mol·K))
    N_Avogadro = 6.02214e23  # Avogadro's number (mol^-1)

    # 2. Pre-calculations for the Mixture
    # Calculate moles and mole fractions
    n_moles_A = N_A / N_Avogadro
    n_moles_B = N_B / N_Avogadro
    n_moles_total = n_moles_A + n_moles_B
    x_A = n_moles_A / n_moles_total
    x_B = n_moles_B / n_moles_total

    # Calculate average molar mass (in kg/mol for SI units)
    M_A_kg = M_A_gmol / 1000.0
    M_B_kg = M_B_gmol / 1000.0
    M_avg = x_A * M_A_kg + x_B * M_B_kg

    # Calculate effective mixture van der Waals parameters using mixing rules
    a_mix = x_A**2 * a_AA + 2 * x_A * x_B * a_AB + x_B**2 * a_BB
    b_mix = x_A * b_AA + x_B * b_BB

    # 3. Define the Ordinary Differential Equation (ODE) for molar density
    def model(z, n_mol):
        """
        Defines the ODE: dn_mol/dz = f(z, n_mol).
        Based on hydrostatic equilibrium and the van der Waals EOS.
        """
        # The term (1 - n_mol * b_mix) must be positive.
        if (1 - n_mol * b_mix) <= 1e-9:
            return np.inf  # Return a large number to signal solver failure
        
        # dP/dn_mol term derived from the van der Waals equation
        dP_dn = (R * T) / (1 - n_mol * b_mix)**2 - 2 * a_mix * n_mol
        
        if abs(dP_dn) < 1e-9:
            return np.inf # Avoid division by zero
        
        # The differential equation for molar density
        dn_mol_dz = -(n_mol * M_avg * g) / dP_dn
        return dn_mol_dz

    # 4. Implement the Shooting Method to find the correct initial density n_mol(0)
    def objective(n0):
        """
        Objective function for the root-finder.
        It calculates the difference between the computed total moles for a given
        initial density n0 and the actual total moles. The root is where this
        difference is zero.
        """
        # Solve the ODE with the guessed initial condition n_mol(0) = n0
        sol = solve_ivp(model, [0, H], [n0], dense_output=True, method='RK45')
        
        # If the solver fails (e.g., n0 is too high), penalize the guess
        if sol.status != 0:
            return n_moles_total 

        # Evaluate the solution on a fine grid for accurate integration
        z_eval = np.linspace(0, H, 200)
        n_mol_profile = sol.sol(z_eval)[0]

        # Ensure the density profile is physically realistic
        if np.any(n_mol_profile < 0):
             return n_moles_total

        # Integrate the profile to find the total calculated moles
        calculated_moles = simps(n_mol_profile, z_eval) * A
        
        return calculated_moles - n_moles_total

    # 5. Solve for n0 and generate the final profile
    # The maximum possible physical density is 1/b_mix
    n_max_physical = 1.0 / b_mix
    
    try:
        # Use a root-finding algorithm to find the correct n_mol(0)
        n0_correct = brentq(objective, 1e-6, n_max_physical * 0.999)
    except ValueError:
        print("Error: Could not determine the initial density n_mol(0).")
        print("The parameters might be physically unrealistic, or the bracketing for the root-finder failed.")
        return

    # Solve the ODE one last time with the correct initial value
    # Evaluate at specific points for output
    z_points_out = np.linspace(0, H, 11)
    final_sol = solve_ivp(model, [0, H], [n0_correct], dense_output=True, t_eval=z_points_out)
    
    z_profile = final_sol.t
    n_mol_profile = final_sol.y[0]
    
    # Calculate the final mass density profile (in kg/m^3)
    rho_profile = n_mol_profile * M_avg

    # 6. Print the results
    print("--- Gas Mixture Properties ---")
    print(f"Average Molar Mass: {M_avg*1000:.4f} g/mol")
    print(f"Mixture vdW 'a_mix': {a_mix:.4f} Pa m^6 mol^-2")
    print(f"Mixture vdW 'b_mix': {b_mix:.4f} m^3 mol^-1")
    print(f"\nSolved Molar Density at z=0: {n0_correct:.4f} mol/m^3\n")
    
    print("--- Final Density Profile ---")
    print("The final relationship is: rho(z) = n_mol(z) * M_avg")
    print(f"where M_avg = {M_avg:.5f} kg/mol, and n_mol(z) is numerically computed.")
    print("-----------------------------------")
    print(" Height (m) | Density (kg/m^3)")
    print("-----------------------------------")
    for z, rho in zip(z_profile, rho_profile):
        print(f"{z:10.2f} | {rho:17.5f}")
    
    # Final answer as requested
    final_answer_val = rho_profile[-1]
    print(f"\n<<<Density at the top (z={H} m) is {final_answer_val:.5f} kg/m^3>>>")

if __name__ == '__main__':
    solve_density_profile()