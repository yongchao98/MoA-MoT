import math

# --- Given Information and Constants ---
# All values are converted to standard SI units (meters, kilograms, seconds, Coulombs, Joules).

# Ionization chamber current in Amperes (C/s)
current = 2.0e-12

# Density of air in kg/m^3
rho_air = 1.293

# Length of ionization chamber in meters
chamber_length = 15.1 / 100

# Photon beam horizontal width at focus in meters
beam_width_h = 0.3 / 1000

# Photon beam vertical width at focus in meters
beam_width_v = 6.0 / 1000

# Total exposure time for a single point on the surface in seconds
exposure_time = 0.02

# Physical constant: Average energy to create an ion pair in air (W_air) per elementary charge (e) in J/C
W_air_over_e = 33.97

# --- Calculation ---

# Step 1: Calculate the cross-sectional area of the beam
beam_area = beam_width_h * beam_width_v

# Step 2: Calculate the mass of the irradiated air in the chamber
mass_air_irradiated = rho_air * beam_area * chamber_length

# Step 3: Calculate the dose rate in Gy/s (J/kg/s)
# Dose Rate = Energy absorbed per second / mass irradiated
# Energy absorbed per second = current (C/s) * W_air/e (J/C)
dose_rate = (current * W_air_over_e) / mass_air_irradiated

# Step 4: Calculate the cumulative dose
cumulative_dose_Gy = dose_rate * exposure_time

# Convert to milliGray (mGy) for readability
cumulative_dose_mGy = cumulative_dose_Gy * 1000

# --- Output the Final Equation and Result ---

print("The cumulative surface dose is calculated as:")
print("Dose = [ (Current * W_air/e) / (Air_Density * Beam_Width * Beam_Height * Chamber_Length) ] * Exposure_Time\n")

print("Final Equation with substituted values:")
# Using scientific notation for clarity where appropriate
print(f"Dose = [ ({current:.1e} C/s * {W_air_over_e} J/C) / ({rho_air} kg/m^3 * {beam_width_h} m * {beam_width_v} m * {chamber_length} m) ] * {exposure_time} s")

# Print the final result
print(f"\nCalculated Cumulative Dose = {cumulative_dose_Gy:.3e} Gy")
print(f"Which is equal to {cumulative_dose_mGy:.4f} mGy")

# The final numerical answer
# print(f"<<<{cumulative_dose_mGy:.4f}>>>")
<<<0.0039>>>