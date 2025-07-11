import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# --- Start of the main script ---

# Constants given in the problem or standard physical values
total_ch4_release_tons = 250000.0  # metric tons
troposphere_mixing_fraction = 0.80
stratosphere_mixing_fraction = 0.20
troposphere_oxidation_y1 = 0.05  # 5% reduction
annual_oxidation_y2_y3 = 0.03    # 3% reduction

# Physical constants
# Using a standard value for the average molecular weight of air.
mw_ch4_g_mol = 16.0
mw_air_g_mol = 28.97
total_atmosphere_mass_kg = 5.1e18

# Step 1: Calculate mass distribution in kg
total_ch4_release_kg = total_ch4_release_tons * 1000.0
mass_in_troposphere_initial = total_ch4_release_kg * troposphere_mixing_fraction
mass_for_stratosphere_total = total_ch4_release_kg * stratosphere_mixing_fraction
stratosphere_annual_addition_kg = mass_for_stratosphere_total / 3.0

# Step 2: Simulate the mass changes year-by-year for 3 years
# At T=0, 80% of methane is in the troposphere. After 1 year, this mass is oxidized.
mass_after_y1_oxidation = mass_in_troposphere_initial * (1 - troposphere_oxidation_y1)
# The first portion of stratospheric methane is then added at the end of Year 1.
mass_end_of_y1 = mass_after_y1_oxidation + stratosphere_annual_addition_kg

# In Year 2, another stratospheric portion is added, and the new total is oxidized.
mass_start_of_y2_process = mass_end_of_y1 + stratosphere_annual_addition_kg
mass_end_of_y2 = mass_start_of_y2_process * (1 - annual_oxidation_y2_y3)

# In Year 3, the final stratospheric portion is added, and the new total is oxidized.
mass_start_of_y3_process = mass_end_of_y2 + stratosphere_annual_addition_kg
final_ch4_mass_kg = mass_start_of_y3_process * (1 - annual_oxidation_y2_y3)

# Step 3: Convert the final mass to an increase in concentration (ppb)
mw_ch4_kg_mol = mw_ch4_g_mol / 1000.0
mw_air_kg_mol = mw_air_g_mol / 1000.0

# Calculate moles of remaining CH4 from the leak
final_ch4_moles = final_ch4_mass_kg / mw_ch4_kg_mol
# Calculate total moles of air in the atmosphere
total_air_moles = total_atmosphere_mass_kg / mw_air_kg_mol

# Calculate the concentration increase in parts-per-billion (ppb)
ppb_increase = (final_ch4_moles / total_air_moles) * 1e9

# --- Output the results and final equation as requested ---
print("Final Calculation Steps:")
print(f"1. Remaining Methane Mass = {final_ch4_mass_kg:.2f} kg")
print(f"2. Final Moles of CH4 from leak = {final_ch4_mass_kg:.2f} kg / {mw_ch4_kg_mol} kg/mol = {final_ch4_moles:.4e} moles")
print(f"3. Total Moles of Air = {total_atmosphere_mass_kg:.1e} kg / {mw_air_kg_mol} kg/mol = {total_air_moles:.4e} moles")
print("-" * 30)
print("Final Equation for ppb Increase:")
# Using the calculated numbers in the final equation string
print(f"ppb increase = ({final_ch4_moles:.4e} moles / {total_air_moles:.4e} moles) * 1,000,000,000")
print("-" * 30)
print(f"The calculated increase in atmospheric methane concentration after 3 years is: {ppb_increase:.4f} ppb")

# --- End of the main script ---

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()

# Print the captured output to the actual console
print(output)
# Print the final answer in the required format
print(f'<<<{ppb_increase:.4f}>>>')