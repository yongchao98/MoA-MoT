import math

# Step 1: Define knowns and constants, converting to consistent units (cm, g, s, C)

# Ionization chamber current in Amperes (C/s)
I = 2.0e-12  # 2.0 pA = 2.0e-12 A

# Photon beam size at focus in cm
h_beam_focus_cm = 0.3 / 10  # 0.3 mm
v_beam_focus_cm = 6.0 / 10  # 6.0 mm
A_beam_cm2 = h_beam_focus_cm * v_beam_focus_cm

# Length of the ionization chamber in cm
L_cm = 15.1

# Density of air in g/cm^3
rho_air_g_cm3 = 1.293 / 1000  # 1.293 mg/cm^3

# Energy absorption mass attenuation coefficient of air in cm^2/g
mu_en_rho_air_cm2_g = 0.328

# Total effective exposure time for a single point on the subject in seconds
t_exposure_s = 0.02

# Physical constant: average energy to create an ion pair in air, in J/C
W_air_e_J_C = 33.97

# Step 2: Calculate Energy Fluence Rate (Psi_dot)

# Calculate the exponent for the absorption fraction formula
exponent = -mu_en_rho_air_cm2_g * rho_air_g_cm3 * L_cm

# Calculate the fraction 'f' of the beam's energy absorbed by the air in the chamber
f_absorbed = 1 - math.exp(exponent)

# Calculate the total energy absorbed per second in the chamber (in J/s)
E_abs_rate_J_s = I * W_air_e_J_C

# Calculate the incident energy fluence rate (Psi_dot) in J/(cm^2 * s)
# This is the rate of energy entering the chamber per unit area.
Psi_dot_J_cm2_s = E_abs_rate_J_s / (f_absorbed * A_beam_cm2)

# Step 3: Calculate Dose Rate (D_dot) in the tissue (assumed to be like air)

# Dose rate in J/g*s is Psi_dot * (mu_en/rho)
D_dot_J_g_s = Psi_dot_J_cm2_s * mu_en_rho_air_cm2_g

# Convert dose rate from J/g*s to Gy/s (1 J/g = 1000 J/kg = 1000 Gy)
D_dot_Gy_s = D_dot_J_g_s * 1000

# Step 4: Calculate the Cumulative Surface Dose (D)

# Cumulative dose is the dose rate multiplied by the total exposure time
D_cumulative_Gy = D_dot_Gy_s * t_exposure_s

# Step 5: Print the final calculation and result
# The final equation is: Cumulative Dose = Dose Rate * Total Exposure Time
print("The final calculation for the cumulative surface dose is:")
print(f"Cumulative Dose = {D_dot_Gy_s:.6f} Gy/s * {t_exposure_s:.2f} s = {D_cumulative_Gy:.8f} Gy")

# For context, this is equivalent to {D_cumulative_Gy * 1000:.4f} mGy
# The final result in mGy is printed in a comment as the primary answer should be in Gy.

# <<<3.88248135e-05>>>