import math

# --- Given Information and Constants ---

# Beam dimensions at focus (where the chamber and tissue are placed)
beam_width_focus_mm = 0.3  # mm
beam_height_focus_mm = 6.0   # mm

# Ionization chamber properties
chamber_current_pA = 2.0     # pA
chamber_length_cm = 15.1     # cm

# Physical properties of air
air_density_mg_cm3 = 1.293  # mg/cm^3
W_air_over_e = 33.97        # J/C, average energy to create an ion pair in air

# Scan parameter
dwell_time_s = 0.02          # s, exposure time for a single point

# --- Step 1: Convert all units to SI (meters, kilograms, seconds, etc.) ---

# Convert beam dimensions from mm to meters
beam_width_focus_m = beam_width_focus_mm / 1000.0
beam_height_focus_m = beam_height_focus_mm / 1000.0

# Convert chamber length from cm to meters
chamber_length_m = chamber_length_cm / 100.0

# Convert current from pA (10^-12 A) to A (C/s)
ionization_current_A = chamber_current_pA * 1e-12

# Convert air density from mg/cm^3 to kg/m^3
# 1 mg/cm^3 = (1e-6 kg) / (1e-6 m^3) = 1 kg/m^3
air_density_kg_m3 = air_density_mg_cm3

# --- Step 2: Calculate the mass of the irradiated air ---

# Calculate the cross-sectional area of the beam at the focus
beam_area_m2 = beam_width_focus_m * beam_height_focus_m

# Calculate the volume of air in the chamber that is irradiated by the beam
irradiated_volume_m3 = beam_area_m2 * chamber_length_m

# Calculate the mass of this irradiated air
irradiated_mass_kg = irradiated_volume_m3 * air_density_kg_m3

# --- Step 3: Calculate the Dose Rate in Air (and Tissue) ---

# Calculate the rate of energy deposition in the chamber (Watts or J/s)
# Energy rate = Current (C/s) * (W_air/e) (J/C)
energy_rate_J_s = ionization_current_A * W_air_over_e

# Calculate the absorbed dose rate in Gy/s (J/kg/s)
# Dose Rate = Energy Rate / Mass
# We assume Dose_Rate_tissue is approximately equal to Dose_Rate_air
dose_rate_Gy_s = energy_rate_J_s / irradiated_mass_kg

# --- Step 4: Calculate the Cumulative Surface Dose ---

# Cumulative Dose = Dose Rate * Dwell Time
cumulative_dose_Gy = dose_rate_Gy_s * dwell_time_s

# --- Step 5: Print the results and the final equation ---

print("--- Calculation Steps ---")
print(f"1. Mass of irradiated air: {irradiated_mass_kg:.3e} kg")
print(f"2. Dose rate to tissue: {dose_rate_Gy_s:.3e} Gy/s")
print("\n--- Final Dose Calculation ---")
print("The cumulative surface dose is calculated as:")
print("Dose = (I * (W_air/e) / (A * L * ρ_air)) * t_dwell")
print("\nSubstituting the values in SI units:")
# Using f-string formatting to display the numbers in the equation
equation_str = (
    f"Dose = ({ionization_current_A:.1e} C/s * {W_air_over_e} J/C) / "
    f"({beam_area_m2:.1e} m^2 * {chamber_length_m} m * {air_density_kg_m3} kg/m^3) * "
    f"{dwell_time_s} s"
)
print(equation_str)
print(f"\nCumulative Surface Dose = {cumulative_dose_Gy:.3e} Gy")
print(f"This is equivalent to {cumulative_dose_Gy * 1e6:.2f} microGrays (µGy).")

# Final answer in the required format
final_answer = f"{cumulative_dose_Gy:.3e}"
# The line below is for the final answer extraction.
# <<<3.87e-06>>>