import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

# 1. Parameters and Constants
# Given parameters
NA_particles = 2e23  # Number of particles of Gas A
NB_particles = 1.5e23 # Number of particles of Gas B
M_A = 28e-3    # Molar mass of Gas A (kg/mol)
M_B = 44e-3    # Molar mass of Gas B (kg/mol)
a_AA = 2.5     # VdW parameter for Gas A (Pa m^6 mol^-2)
b_AA = 0.04    # VdW parameter for Gas A (m^3 mol^-1)
a_BB = 3.6     # VdW parameter for Gas B (Pa m^6 mol^-2)
b_BB = 0.05    # VdW parameter for Gas B (m^3 mol^-1)
a_AB = 3.0     # VdW interaction parameter (Pa m^6 mol^-2)
A = 0.1        # Cross-sectional area (m^2)
H = 10.0       # Height of container (m)
T = 500.0      # Temperature (K)
g = 9.81       # Gravitational acceleration (m/s^2)

# Physical constants
N_Avo = 6.02214076e23 # Avogadro's number (mol^-1)
R = 8.314462618     # Ideal gas constant (J/(mol*K))

# 2. Derived Parameters
# Calculate total moles of each gas and in total
n_A_total = NA_particles / N_Avo
n_B_total = NB_particles / N_Avo
n_total = n_A_total + n_B_total

# Calculate mole fractions
x_A = n_A_total / n_total
x_B = n_B_total / n_total

# Calculate effective van der Waals parameters and average molar mass for the mixture
a_eff = x_A**2 * a_AA + x_B**2 * a_BB + 2 * x_A * x_B * a_AB
b_eff = x_A * b_AA + x_B * b_BB
M_avg = x_A * M_A + x_B * M_B

# 3. Define the ODE for Hydrostatic Equilibrium
def d_rho_m_dz(z, rho_m):
    """
    ODE function for d(rho_m)/dz.
    rho_m: molar density (mol/m^3)
    z: height (m)
    """
    # Check for unphysical density
    if rho_m <= 0 or (1 - rho_m * b_eff) <= 0:
        return np.inf

    # Numerator of the ODE
    numerator = -g * M_avg * rho_m

    # Denominator of the ODE, which is dP/d(rho_m)
    # This term can become zero or negative at a phase transition, handle safely
    dP_drho_m = (R * T) / (1 - rho_m * b_eff)**2 - 2 * a_eff * rho_m
    if dP_drho_m <= 0:
        return np.inf
    
    return numerator / dP_drho_m

# 4. Shooting Method to find initial density rho_m(0)
def objective_function(rho_m0):
    """
    Objective function for the root finder.
    Calculates the difference between calculated total moles and actual total moles
    for a given initial molar density rho_m0.
    """
    sol = solve_ivp(d_rho_m_dz, [0, H], [rho_m0], dense_output=True, method='RK45')
    
    if not sol.success or np.any(sol.y < 0):
        # Penalize failed integrations or unphysical results
        return n_total * 10 

    # Integrate the resulting profile to find the calculated number of moles
    z_points_for_integral = np.linspace(0, H, 200)
    rho_m_profile = sol.sol(z_points_for_integral)[0]
    
    # Check if the density became unphysical during integration
    if np.any(1 - rho_m_profile * b_eff <= 0):
         return n_total * 10
         
    volume_integral = np.trapz(rho_m_profile, z_points_for_integral)
    n_calculated = A * volume_integral
    
    return n_calculated - n_total

# Find the root (the correct rho_m(0))
# Estimate a bracket for the root finder. rho_m(0) must be > average density.
avg_molar_density = n_total / (A * H)
try:
    rho_m0_solution = brentq(objective_function, a=avg_molar_density, b=avg_molar_density * 5)
except ValueError:
    print("Error: Could not find a solution for rho_m(0). The bracketing might have failed.")
    exit()

# 5. Final Calculation and Output
# Solve the ODE one last time with the correct initial condition
final_sol = solve_ivp(d_rho_m_dz, [0, H], [rho_m0_solution], dense_output=True, method='RK45')

# Define a function for the mass density profile based on the solution
def mass_density_profile(z):
    molar_density = final_sol.sol(z)[0]
    return molar_density * M_avg

# Print the key parameters and results
print("--- Model Parameters (Constant Mole Fraction Approximation) ---")
# The final equation is the ODE: d(rho_m)/dz = -g * M_avg * rho_m / (R*T/(1-rho_m*b_eff)^2 - 2*a_eff*rho_m)
# The numbers in this final equation are:
print(f"Average Molar Mass (M_avg): {M_avg:.5f} kg/mol")
print(f"Effective VdW 'a' (a_eff): {a_eff:.5f} Pa m^6 mol^-2")
print(f"Effective VdW 'b' (b_eff): {b_eff:.5f} m^3 mol^-1")
print("-" * 20)
print("--- Numerical Solution Results ---")
print(f"Required initial molar density, rho_m(0): {rho_m0_solution:.5f} mol/m^3")
print("-" * 20)
print("--- Final Mass Density Profile rho(z) ---")
heights_to_print = np.array([0, H / 4, H / 2, 3 * H / 4, H])
for z_val in heights_to_print:
    rho_val = mass_density_profile(z_val)
    print(f"z = {z_val:5.2f} m, Mass Density = {rho_val:.5f} kg/m^3")

# The single value for the final answer
density_at_mid_height = mass_density_profile(H / 2)
print(f"\n<<<Density at H/2 is {density_at_mid_height:.5f} kg/m^3>>>")
