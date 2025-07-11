import math

# --- Setup: Constants and Parameters in SI units ---

# Physical constants
K_B = 1.380649e-23  # Boltzmann constant in J/K
E_CHARGE = 1.6021766e-19   # Elementary charge in Coulombs

# Given parameters
T_KELVIN = 3000.0  # Temperature in Kelvin
P_TORR = 10.0      # Pressure in Torr
I_W_CM2 = 10.0     # Intensity in W/cm^2

# Derived and literature values
E_ION_H_EV = 13.6  # Ionization energy of hydrogen in eV
SIGMA_THRESHOLD = 6.3e-22 # Photoionization cross-section at 13.6 eV in m^2

# --- Unit Conversions to SI ---
P_PASCAL = P_TORR * 133.322   # Pressure from Torr to Pa
I_W_M2 = I_W_CM2 * 1e4        # Intensity from W/cm^2 to W/m^2
E_ION_H_JOULES = E_ION_H_EV * E_CHARGE # Ionization energy from eV to Joules

# --- Step-by-Step Calculation ---

# 1. Calculate density of hydrogen atoms (n_H) using the Ideal Gas Law
n_H = P_PASCAL / (K_B * T_KELVIN)
print(f"1. Calculating the density of hydrogen atoms (n_H)...")
print(f"   Using the Ideal Gas Law: n_H = P / (k_B * T)")
print(f"   -> n_H = {P_PASCAL:.2f} Pa / ({K_B:.4e} J/K * {T_KELVIN:.0f} K)")
print(f"   n_H = {n_H:.3e} m^-3\n")

# 2. Calculate the energy of an incident photon (E_photon)
# The radiation frequency ω ~ e^2/(ħ*a_B) corresponds to an energy of 2 * Ry = 27.2 eV
E_photon = 2 * E_ION_H_JOULES
print(f"2. Calculating the energy per photon (E_photon)...")
print(f"   The photon energy is given as twice the hydrogen ionization energy.")
print(f"   -> E_photon = 2 * {E_ION_H_JOULES:.3e} J")
print(f"   E_photon = {E_photon:.3e} J\n")

# 3. Calculate the photon flux (Phi)
Phi = I_W_M2 / E_photon
print(f"3. Calculating the photon flux (Phi)...")
print(f"   Flux is Intensity divided by Energy per photon: Phi = I / E_photon")
print(f"   -> Phi = {I_W_M2:.1e} W/m^2 / {E_photon:.3e} J")
print(f"   Phi = {Phi:.3e} m^-2 s^-1\n")

# 4. Estimate the photoionization cross-section (sigma_ph)
# The cross-section decreases with energy: σ(E) ≈ σ_threshold * (E_ion / E)^3
E_ratio = E_ION_H_JOULES / E_photon  # This will be 0.5
sigma_ph = SIGMA_THRESHOLD * (E_ratio**3)
print(f"4. Estimating the photoionization cross-section (sigma_ph)...")
print(f"   Using the approximation: sigma_ph ≈ sigma_threshold * (E_ion / E_photon)^3")
print(f"   -> sigma_ph ≈ {SIGMA_THRESHOLD:.1e} m^2 * ({E_ratio})**3")
print(f"   sigma_ph = {sigma_ph:.3e} m^2\n")

# 5. Estimate the recombination coefficient (alpha)
# Using a standard formula for radiative recombination in a plasma
# alpha ≈ 2.6e-13 * (T/10^4 K)^-0.7 cm^3/s
alpha_cm3_s = 2.6e-13 * (T_KELVIN / 1e4)**(-0.7)
alpha_m3_s = alpha_cm3_s * 1e-6  # Convert cm^3/s to m^3/s
print(f"5. Estimating the recombination coefficient (alpha)...")
print(f"   Using an empirical formula for T = {T_KELVIN:.0f} K.")
print(f"   alpha ≈ {alpha_m3_s:.3e} m^3 s^-1\n")

# 6. Solve for electron density (n_e) in steady state
# Rate_ionization = Rate_recombination => n_H * sigma_ph * Phi = alpha * n_e^2
numerator = n_H * sigma_ph * Phi
denominator = alpha_m3_s
n_e = math.sqrt(numerator / denominator)

print(f"6. Solving for the steady-state photoelectron density (n_e)...")
print(f"   The final equation is: n_e = sqrt( (n_H * sigma_ph * Phi) / alpha )")
print(f"   The values for the equation are:")
print(f"   n_H      = {n_H:.3e} m^-3")
print(f"   sigma_ph = {sigma_ph:.3e} m^2")
print(f"   Phi      = {Phi:.3e} m^-2 s^-1")
print(f"   alpha    = {alpha_m3_s:.3e} m^3 s^-1")
print(f"\n   n_e = sqrt( ({n_H:.3e} * {sigma_ph:.3e} * {Phi:.3e}) / {alpha_m3_s:.3e} )")
print(f"\nEstimated density of photoelectrons:")
print(f"n_e ≈ {n_e:.3e} m^-3")
<<<3.129e+20>>>