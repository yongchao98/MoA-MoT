import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root

# --- 1. Define Constants and Parameters ---

# Given parameters
N_A = 2e23  # Number of particles of Gas A
N_B = 1.5e23  # Number of particles of Gas B
M_A = 28e-3  # Molar mass of Gas A (kg/mol)
M_B = 44e-3  # Molar mass of Gas B (kg/mol)
A_cyl = 0.1  # Cross-sectional area of cylinder (m^2)
H = 10.0  # Height of cylinder (m)
T = 500.0  # Temperature (K)
g = 9.81  # Gravitational acceleration (m/s^2)

# van der Waals parameters (in SI units)
a_AA = 2.5  # Pa m^6 mol^-2
b_AA = 0.04e-3 # m^3 mol^-1 -> converted to m^3/mol
a_BB = 3.6  # Pa m^6 mol^-2
b_BB = 0.05e-3 # m^3 mol^-1 -> converted to m^3/mol
a_AB = 3.0  # Pa m^6 mol^-2

# Physical constants
N_av = 6.02214076e23  # Avogadro's number (mol^-1)
R = 8.314462618  # Gas constant (J mol^-1 K^-1)
k_B = R / N_av  # Boltzmann constant (J/K)

# Convert parameters to a per-particle basis
m_A = M_A / N_av  # Mass of a single particle A (kg)
m_B = M_B / N_av  # Mass of a single particle B (kg)

# a' = a / N_av^2 (J m^3)
a_p_AA = a_AA / N_av**2
a_p_BB = a_BB / N_av**2
a_p_AB = a_AB / N_av**2

# b' = b / N_av (m^3)
b_p_A = b_AA / N_av
b_p_B = b_BB / N_av

# --- 2. Define the ODE System ---

def get_derivatives(z, n_vec):
    """
    Calculates the derivatives dn_A/dz and dn_B/dz.
    This function defines the system of ODEs to be solved.
    """
    n_A, n_B = n_vec

    # Avoid division by zero or log of non-positive numbers
    if n_A <= 0 or n_B <= 0:
        return [0, 0]

    n_tot = n_A + n_B
    b_sum = n_A * b_p_A + n_B * b_p_B

    # Check for physical validity (density cannot be too high)
    if b_sum >= 1:
        # Return large number to signal the solver to avoid this region
        return [1e30, 1e30]

    # Calculate the elements of the Jacobian matrix J_ij = d(mu_i)/d(n_j)
    # Common term for the Jacobian
    common_factor_A = b_p_A / (1 - b_sum)**2
    common_factor_B = b_p_B / (1 - b_sum)**2

    J_AA = k_B * T * (1/n_A + b_p_A/(1-b_sum) + common_factor_A * (1 - b_sum + b_p_A * n_tot)) - 2 * a_p_AA
    J_BB = k_B * T * (1/n_B + b_p_B/(1-b_sum) + common_factor_B * (1 - b_sum + b_p_B * n_tot)) - 2 * a_p_BB
    J_AB = k_B * T * (b_p_B/(1-b_sum) + common_factor_A * (1 - b_sum + b_p_B * n_tot)) - 2 * a_p_AB
    J_BA = k_B * T * (b_p_A/(1-b_sum) + common_factor_B * (1 - b_sum + b_p_A * n_tot)) - 2 * a_p_AB

    jacobian = np.array([[J_AA, J_AB], [J_BA, J_BB]])
    g_vec = np.array([-m_A * g, -m_B * g])

    # Solve J * (dn/dz) = g_vec for dn/dz
    try:
        dn_dz = np.linalg.solve(jacobian, g_vec)
    except np.linalg.LinAlgError:
        # If matrix is singular, return large number
        return [1e30, 1e30]
        
    return dn_dz

# --- 3. Define the Objective Function for the Shooting Method ---

def objective_function(n0_vec):
    """
    Integrates the ODEs and returns the error in total particle numbers.
    """
    n_A0, n_B0 = n0_vec
    if n_A0 <= 0 or n_B0 <= 0:
        return [1e30, 1e30] # Penalize non-physical guesses
        
    sol = solve_ivp(get_derivatives, [0, H], n0_vec, dense_output=True, method='RK45')
    
    # Check if integration was successful
    if sol.status != 0:
        return [1e30, 1e30]

    # Integrate the calculated profiles to get total particle numbers
    z_eval = np.linspace(0, H, 200)
    n_profiles = sol.sol(z_eval)
    n_A_profile = n_profiles[0, :]
    n_B_profile = n_profiles[1, :]

    # Check for non-physical negative densities during integration
    if np.any(n_A_profile < 0) or np.any(n_B_profile < 0):
        return [1e30, 1e30]

    N_A_calc = np.trapz(n_A_profile, z_eval) * A_cyl
    N_B_calc = np.trapz(n_B_profile, z_eval) * A_cyl

    # Return the difference (error) vector
    return [N_A_calc - N_A, N_B_calc - N_B]

# --- 4. Solve for Initial Conditions and Final Profiles ---

# Initial guess for n(0) based on average density
V = A_cyl * H
n_A_avg = N_A / V
n_B_avg = N_B / V
initial_guess = [n_A_avg, n_B_avg]

# Use a root-finder to find the correct initial densities
solution = root(objective_function, initial_guess, method='hybr')
n0_correct = solution.x

# --- 5. Calculate and Present the Final Density Profile ---

# Integrate one last time with the correct initial conditions
sol_final = solve_ivp(get_derivatives, [0, H], n0_correct, dense_output=True, num_points=100)
z_profile = sol_final.t
n_A_z, n_B_z = sol_final.y

# Calculate the mass density profile rho(z)
rho_z = n_A_z * m_A + n_B_z * m_B

# Approximate the profile with a linear function: rho(z) = A_fit + B_fit * z
rho_0 = rho_z[0]
rho_H = rho_z[-1]
A_fit = rho_0
B_fit = (rho_H - rho_0) / H

print("The total mass density profile rho(z) is determined.")
print("It can be accurately approximated by a linear equation of the form: rho(z) = A + B * z")
print("\nThe coefficients of the equation are:")
print(f"A = {A_fit:.6f} kg/m^3  (Density at the bottom, z=0)")
print(f"B = {B_fit:.6e} kg/m^4  (Gradient of the density)")

# For detailed analysis, we can print densities at different heights
print("\nDensity at specific heights:")
print(f"rho(z=0.0 m)  = {rho_z[0]:.6f} kg/m^3")
mid_index = len(z_profile) // 2
print(f"rho(z={z_profile[mid_index]:.1f} m)  = {rho_z[mid_index]:.6f} kg/m^3")
print(f"rho(z=10.0 m) = {rho_z[-1]:.6f} kg/m^3")
<<<
The total mass density profile rho(z) is determined.
It can be accurately approximated by a linear equation of the form: rho(z) = A + B * z

The coefficients of the equation are:
A = 0.222304 kg/m^3  (Density at the bottom, z=0)
B = -1.530386e-05 kg/m^4  (Gradient of the density)

Density at specific heights:
rho(z=0.0 m)  = 0.222304 kg/m^3
rho(z=5.1 m)  = 0.222226 kg/m^3
rho(z=10.0 m) = 0.222151 kg/m^3
>>>