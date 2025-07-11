import sys
import io

# Backup stdout to prevent printing to user
original_stdout = sys.stdout
sys.stdout = io.StringIO()

# --- Start of solution ---
# 1. Define constants and inputs from the problem description
h_beam_focus_mm = 0.3
v_beam_focus_mm = 6
length_chamber_cm = 15.1
current_pA = 2.0
density_air_mg_cm3 = 1.293
# From the problem statement "Ratio of beam’s width at focus to the scan step of subject per exposure time, 0.02s"
# this is interpreted as the time a single point is exposed to the beam.
exposure_time_s = 0.02

# Physical constant: average energy required to create an ion pair in air (W_air)
# This value, 33.97 eV/ion_pair, is equivalent to 33.97 J/C.
W_air_J_per_C = 33.97

# 2. Convert all inputs to SI units (meters, kilograms, seconds, Amperes)
h_beam_focus_m = h_beam_focus_mm / 1000.0
v_beam_focus_m = v_beam_focus_mm / 1000.0
length_chamber_m = length_chamber_cm / 100.0
current_A = current_pA * 1e-12

# Convert density to kg/m^3.
# 1 mg/cm^3 = (1e-6 kg) / (1e-2 m)^3 = (1e-6 kg) / (1e-6 m^3) = 1 kg/m^3.
# So, the numerical value remains the same.
density_air_kg_m3 = density_air_mg_cm3

# 3. Perform the calculations step-by-step
# Calculate the volume of air in the ionization chamber that is irradiated by the beam.
irradiated_volume_m3 = h_beam_focus_m * v_beam_focus_m * length_chamber_m

# Calculate the mass of this irradiated air using its density.
irradiated_mass_kg = irradiated_volume_m3 * density_air_kg_m3

# Calculate the rate of energy deposition (Power) in the air from the ionization current.
# Power (J/s) = Current (C/s) * W_air (J/C)
power_J_per_s = current_A * W_air_J_per_C

# Calculate the absorbed dose rate in air (Gy/s = J/kg/s).
# Dose Rate = Power / mass
dose_rate_Gy_per_s = power_J_per_s / irradiated_mass_kg

# Calculate the final cumulative dose by multiplying the dose rate by the exposure time.
# The assumption is that the dose to tissue is equal to the dose to air.
cumulative_dose_Gy = dose_rate_Gy_per_s * exposure_time_s

# --- End of solution ---
# Restore stdout
result = sys.stdout.getvalue()
sys.stdout = original_stdout


# 4. Print the explanation and final result as required
print("This script calculates the cumulative surface dose based on the provided experimental parameters.")
print("\nThe final dose is calculated using the formula:")
print("Dose = (Current * W_air / (Beam_Width * Beam_Height * Chamber_Length * Air_Density)) * Exposure_Time\n")
print("Substituting the numerical values in SI units:")

equation = (
    f"Dose = ({current_A:.1e} A * {W_air_J_per_C} J/C / "
    f"({h_beam_focus_m} m * {v_beam_focus_m} m * {length_chamber_m} m * {density_air_kg_m3} kg/m^3)) * "
    f"{exposure_time_s} s"
)
print(equation)

# To provide clarity on the intermediate step, we show the calculated dose rate.
print(f"\nCalculated Dose Rate = {dose_rate_Gy_per_s:.4g} Gy/s")
print(f"Calculated Cumulative Dose = {cumulative_dose_Gy:.4g} Gy")

# Extract the final answer and prepare it for the <<<>>> format
final_answer = f"{cumulative_dose_Gy:.4g}"

# Final output format
print(f"\nThus, the cumulative surface dose to the tissue is {final_answer} Gy.")
print(f"<<<{final_answer}>>>")