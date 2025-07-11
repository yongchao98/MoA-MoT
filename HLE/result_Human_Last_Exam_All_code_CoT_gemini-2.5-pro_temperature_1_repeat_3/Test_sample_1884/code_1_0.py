import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# --- Step 1: Define constants and initial values ---
total_ch4_release_tons = 250000
ch4_molar_mass_g_mol = 16.0
atmosphere_total_mass_kg = 5.1e18
# Average molar mass of air in g/mol. Standard value is ~28.97 g/mol.
air_molar_mass_g_mol = 29.0

# Convert units for consistency (to kg and kg/mol)
total_ch4_release_kg = total_ch4_release_tons * 1000  # tons to kg
ch4_molar_mass_kg_mol = ch4_molar_mass_g_mol / 1000  # g/mol to kg/mol
air_molar_mass_kg_mol = air_molar_mass_g_mol / 1000  # g/mol to kg/mol

# --- Step 2: Track mass reduction over 3 years ---
print("--- Methane Mass Balance Calculation ---")

# Initial mass
mass_y0 = total_ch4_release_kg
print(f"Initial excess CH4 mass: {mass_y0:e} kg")

# Year 1
tropospheric_mix_fraction_y1 = 0.80
oxidation_rate_y1 = 0.05
mass_in_troposphere_y1 = mass_y0 * tropospheric_mix_fraction_y1
mass_reduction_y1 = mass_in_troposphere_y1 * oxidation_rate_y1
mass_y1 = mass_y0 - mass_reduction_y1
print(f"Mass reduction in Year 1 (5% of 80% portion): {mass_reduction_y1:e} kg")
print(f"Mass at end of Year 1: {mass_y1:e} kg")

# Year 2
oxidation_rate_subsequent = 0.03
mass_reduction_y2 = mass_y1 * oxidation_rate_subsequent
mass_y2 = mass_y1 - mass_reduction_y2
print(f"Mass reduction in Year 2 (3% of remainder): {mass_reduction_y2:e} kg")
print(f"Mass at end of Year 2: {mass_y2:e} kg")

# Year 3
mass_reduction_y3 = mass_y2 * oxidation_rate_subsequent
final_mass_ch4_kg = mass_y2 - mass_reduction_y3
print(f"Mass reduction in Year 3 (3% of remainder): {mass_reduction_y3:e} kg")
print(f"Final excess CH4 mass after 3 years: {final_mass_ch4_kg:e} kg")
print("-" * 36)

# --- Step 3 & 4: Convert to moles ---
print("\n--- Molar Calculation ---")
final_moles_ch4 = final_mass_ch4_kg / ch4_molar_mass_kg_mol
print(f"Final excess CH4 moles: {final_moles_ch4:e} mol")

total_moles_air = atmosphere_total_mass_kg / air_molar_mass_kg_mol
print(f"Total moles of air in atmosphere: {total_moles_air:e} mol")
print("-" * 25)

# --- Step 5: Calculate final concentration increase in ppb ---
print("\n--- Final Concentration Calculation ---")
ppb_increase = (final_moles_ch4 / total_moles_air) * 1e9

# Print the final equation with all numbers
print("Increase in ppb = (Final Moles CH4 / Total Moles Air) * 1,000,000,000")
print(f"Increase in ppb = ({final_moles_ch4:.5e} / {total_moles_air:.5e}) * 1,000,000,000")
print(f"Result: {ppb_increase:.2f} ppb")

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()
# Print the output to the actual console
print(output)
print(f'<<<{ppb_increase:.2f}>>>')