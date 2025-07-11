import math

# --- 1. Define Constants and Given Parameters (in SI units) ---

# Physical Constants
k_B = 1.380649e-23  # Boltzmann constant in J/K
Ry_eV = 13.6057     # Rydberg energy in eV
e_charge = 1.60218e-19 # Elementary charge in C
Ry_J = Ry_eV * e_charge # Rydberg energy in Joules

# Given Parameters
P_torr = 10.0
T_K = 3000.0
I_W_cm2 = 10.0

# Conversions to SI units
P_Pa = P_torr * 133.322         # Pressure in Pascals
I_W_m2 = I_W_cm2 * (100**2)     # Intensity in W/m^2

# --- 2. Calculate Intermediate Values ---

# Calculate initial hydrogen atom density using the ideal gas law: n_H = P / (k_B * T)
n_H = P_Pa / (k_B * T_K)

# Calculate photon energy: E_photon = 2 * Ry
E_photon = 2 * Ry_J

# Calculate photon flux: Phi = I / E_photon
Phi = I_W_m2 / E_photon

# Estimate photoionization cross-section at E = 2*Ry.
# sigma_pi(E) approx sigma_0 * (Ry/E)^3, where sigma_0 is cross-section at threshold (E=Ry)
# sigma_0 for Hydrogen is approx 6.3e-18 cm^2 = 6.3e-22 m^2
sigma_0 = 6.3e-22 # m^2
sigma_pi = sigma_0 * (Ry_J / E_photon)**3

# Estimate radiative recombination coefficient alpha_r for T_e = 3000 K.
# alpha_r(T_e) approx 2.7e-13 * (T_e[K] / 10^4)^-0.7 cm^3/s
# We assume electron temperature T_e is equal to the gas temperature T_K.
alpha_r_cm3_s = 2.7e-13 * (T_K / 1e4)**(-0.7)
alpha_r_m3_s = alpha_r_cm3_s * 1e-6 # Convert from cm^3/s to m^3/s

# --- 3. Calculate Final Electron Density ---

# In steady state, production rate = loss rate: n_H * sigma_pi * Phi = alpha_r * n_e^2
# Solve for n_e: n_e = sqrt((n_H * sigma_pi * Phi) / alpha_r)
n_e_squared = (n_H * sigma_pi * Phi) / alpha_r_m3_s
n_e = math.sqrt(n_e_squared)

# --- 4. Output the result ---
# The final equation is n_e = sqrt((n_H * sigma_pi * Phi) / alpha_r)
# The following print statement shows the numbers used in this equation and the final result.
print("Final Equation with calculated values:")
print(f"n_e = sqrt( ( {n_H:.3e} [m^-3] * {sigma_pi:.3e} [m^2] * {Phi:.3e} [m^-2 s^-1] ) / {alpha_r_m3_s:.3e} [m^3 s^-1] )")
print("\nResult:")
print(f"The estimated density of photoelectrons is {n_e:.3e} m^-3.")