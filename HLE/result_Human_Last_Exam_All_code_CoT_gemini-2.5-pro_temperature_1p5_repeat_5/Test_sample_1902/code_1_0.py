import math

# Step 1: Define given parameters and constants.
# Ionization chamber current in Amperes (C/s)
i_chamber_A = 2.0e-12
# Photon beam width at focus in cm
beam_w_cm = 0.3 / 10.0
# Photon beam height at focus in cm
beam_h_cm = 6.0 / 10.0
# Length of ionization chamber in cm
L_chamber_cm = 15.1
# Density of air in g/cm^3
rho_air_g_cm3 = 1.293 / 1000.0
# Effective exposure time for a point on the subject in seconds
t_exposure_s = 0.02
# Average energy to produce an ion pair in air (W_air/e) in J/C
W_air_per_e_J_C = 33.97

# Step 2: Calculate the mass of the air being ionized in grams.
# Beam area in cm^2
beam_area_cm2 = beam_w_cm * beam_h_cm
# Volume of air in the beam's path within the chamber in cm^3
volume_air_cm3 = beam_area_cm2 * L_chamber_cm
# Mass of the ionized air in grams
m_air_g = volume_air_cm3 * rho_air_g_cm3
# Convert mass to kilograms for dose calculation
m_air_kg = m_air_g / 1000.0

# Step 3: Calculate the rate of energy absorption (Power) in J/s.
power_J_s = i_chamber_A * W_air_per_e_J_C

# Step 4: Determine the dose rate in Gray per second (Gy/s).
# 1 Gy = 1 J/kg
dose_rate_Gy_s = power_J_s / m_air_kg

# Step 5: Calculate the final cumulative dose in Gray (Gy).
cumulative_dose_Gy = dose_rate_Gy_s * t_exposure_s

# Step 6: Print the final equation with all numbers and the result.
# We will show the full calculation for the final dose.
# Cumulative Dose = (Power / Mass) * Exposure Time
# Cumulative Dose = ((Current * W/e) / (Width * Height * Length * Density)) * Exposure Time
print("The cumulative surface dose calculation is based on the formula:")
print("Dose = ((Current * (W/e)) / (Beam_Width * Beam_Height * Chamber_Length * Air_Density)) * Exposure_Time\n")

print("Substituting the given values into the equation:")
# We use scientific notation for clarity in the equation string.
print(f"Dose = (({i_chamber_A:.1e} C/s * {W_air_per_e_J_C} J/C) / ({beam_w_cm} cm * {beam_h_cm} cm * {L_chamber_cm} cm * {rho_air_g_cm3:.4e} g/cm^3 / 1000 g/kg)) * {t_exposure_s} s")

# Print the intermediate steps
print(f"\nEnergy Absorption Rate = {i_chamber_A:.1e} C/s * {W_air_per_e_J_C} J/C = {power_J_s:.3e} J/s")
print(f"Mass of Ionized Air = {beam_area_cm2:.3f} cm^2 * {L_chamber_cm} cm * {rho_air_g_cm3:.4e} g/cm^3 = {m_air_g:.3e} g = {m_air_kg:.3e} kg")
print(f"Dose Rate = {power_J_s:.3e} J/s / {m_air_kg:.3e} kg = {dose_rate_Gy_s:.5f} Gy/s")

# Print the final result
print(f"\nCumulative Dose = {dose_rate_Gy_s:.5f} Gy/s * {t_exposure_s} s = {cumulative_dose_Gy:.3e} Gy")
print(f"The cumulative surface dose to the tissue is {cumulative_dose_Gy:.3e} Gray.")

<<<3.866e-06>>>